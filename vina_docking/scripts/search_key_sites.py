import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from pathlib import Path
from glob import glob
from tqdm import tqdm
from typing import Optional
import subprocess
import warnings
import requests
from lxml import etree
import multitasking
import numpy as np
from itertools import chain
from retry import retry
import os

warnings.filterwarnings("ignore")


DB_PATH = "/home/guolj/ws/database/uniprot_sprot/sprot"
MMSEQS = "/home/guolj/miniconda3/envs/docking/bin/mmseqs"
SEARCH_FORMAT = [
    "query",
    "target",
    "evalue",
    "pident",
    "qcov",
    "tcov",
    "qseq",
    "tseq",
    "qaln",
    "taln",
    "qstart",
    "qend",
    "tstart",
    "tend",
]
XML_DB = "/home/guolj/ws/database/uniprot_xmls"
NAMESPACE = {"ns": "http://uniprot.org/uniprot"}


@retry(tries=3, delay=1)
@multitasking.task
def download_xml(uniprot_id: str):
    if Path(XML_DB, f"{uniprot_id}.xml").exists():
        return
    ret = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml").text
    Path(XML_DB, f"{uniprot_id}.xml").write_text(ret)


def get_features(uniprot_id: str):
    xml_path = Path(XML_DB, f"{uniprot_id}.xml")
    assert xml_path.exists(), f"{uniprot_id} not exists"
    ret = xml_path.read_text()
    xml_tree = etree.fromstring(ret.encode("utf-8"))
    features = xml_tree.findall(".//ns:feature", NAMESPACE)
    feature_positions = []
    for feature in features:
        feature_type = feature.attrib.get("type")
        feature_id = feature.attrib.get("id")

        # 获取description
        description_element = feature.find(".//ns:description", NAMESPACE)
        feature_description = (
            description_element.text if description_element is not None else ""
        )

        # 获取position
        position_element = feature.find(".//ns:location/ns:begin", NAMESPACE)
        if position_element is not None:
            feature_position = position_element.attrib.get("position", "unknown")
            end_element = feature.find(".//ns:location/ns:end", NAMESPACE)
            if end_element is not None:
                feature_position = f"{feature_position}-{end_element.attrib.get('position', 'unknown')}"
        else:
            position_element = feature.find(".//ns:location/ns:position", NAMESPACE)
            feature_position = (
                position_element.attrib.get("position")
                if position_element is not None
                else ""
            )
        feature_positions.append((feature_type, feature_position))
    return feature_positions


def agg_key_sites(in_series: pd.Series):
    return ";".join([i for i in in_series.values if i])


def search_key_sites(docking_results_path: str, outdir: Optional[str] = None, n_workers: int = os.cpu_count(), multitasking_engine: str = "process"):
    multitasking.set_max_threads(n_workers)
    multitasking.set_engine(multitasking_engine)
    docking_results_dir = Path(docking_results_path).parent
    docking_results = pd.read_csv(docking_results_path)
    all_receptors = glob(str(Path(docking_results_dir, "receptors_pdbqt/*.pdbqt")))

    pdb_parser = PDBParser()
    receptor_id_to_sequence = {
        Path(receptor).stem: {
            chain_id: "".join(
                seq1(residue.get_resname())
                for residue in residues
                if residue.get_resname() != "HOH"
            )
            for chain_id, residues in pdb_parser.get_structure("", receptor)[
                0
            ].child_dict.items()
        }
        for receptor in tqdm(all_receptors, desc="get receptor sequences")
    }

    search_results_dir = (
        str(Path(docking_results_dir, "search_results"))
        if outdir is None
        else str(Path(outdir))
    )
    Path(search_results_dir).mkdir(exist_ok=True, parents=True)

    query_path = str(Path(search_results_dir, "seqs.fasta"))
    with open(query_path, "w") as f:
        for receptor_id, chains in receptor_id_to_sequence.items():
            for chain_id, sequence in chains.items():
                f.write(f">{receptor_id}::{chain_id}\n{sequence}\n")
    mmseqs_ret = subprocess.run(
        " ".join(
            [
                MMSEQS,
                "easy-search",
                query_path,
                DB_PATH,
                str(Path(search_results_dir, "alnRes.m8")),
                str(Path(search_results_dir, "tmp")),
                "--format-output",
                f"\"{','.join(SEARCH_FORMAT)}\"",
            ]
        ),
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    search_result = pd.read_csv(
        str(Path(search_results_dir, "alnRes.m8")),
        sep="\t",
        names=SEARCH_FORMAT,
    ).sort_values(
        ["query", "evalue", "pident", "qcov", "tcov"],
        ascending=[True, True, False, False, False],
    )

    for uniprot_id in tqdm(
        search_result["target"].unique(), desc="download uniprot xml"
    ):
        download_xml(uniprot_id)

    search_result_features = {}
    search_result_key_sites = {}
    key_sites_to_feature = {}
    for uniprot_id in tqdm(search_result["target"].unique(), desc="get features"):
        try:
            features = get_features(uniprot_id)
            search_result_features[uniprot_id] = " / ".join(
                [": ".join(i) for i in features]
            )
            search_result_key_sites[uniprot_id] = []
            if uniprot_id not in key_sites_to_feature:
                key_sites_to_feature[uniprot_id] = {}
            for curr_feat, curr_sites in features:
                if "site" in curr_feat:
                    if "-" not in curr_sites:
                        search_result_key_sites[uniprot_id].append(int(curr_sites))
                        key_sites_to_feature[uniprot_id][int(curr_sites)] = curr_feat
                    else:
                        start, end = list(map(int, curr_sites.split("-")))
                        search_result_key_sites[uniprot_id].extend(
                            list(range(start, end + 1))
                        )
                        for site in range(start, end + 1):
                            key_sites_to_feature[uniprot_id][site] = curr_feat
        except Exception as e:
            print(e)
            print(uniprot_id)
            print(str(Path(XML_DB, f"{uniprot_id}.xml")))
            raise Exception

    search_result["target_features"] = search_result["target"].map(
        search_result_features
    )
    search_result["target_key_sites"] = search_result["target"].map(
        lambda x: ";".join(sorted(set(map(str, search_result_key_sites[x])), key=lambda x: int(x)))
    )

    query_key_sites = []
    query_key_sites_features = []
    for tid, qseq, qaln_seq, qstart, tseq, taln_seq, tstart, t_key_sites in tqdm(
        search_result[
            [
                "target",
                "qseq",
                "qaln",
                "qstart",
                "tseq",
                "taln",
                "tstart",
                "target_key_sites",
            ]
        ].values,
        desc="get query key sites",
    ):
        curr_tid_to_features = key_sites_to_feature[tid]
        curr_query_key_sites = []
        if t_key_sites:
            t_key_sites = list(map(int, t_key_sites.split(";")))
            t_key_sites_max = max(t_key_sites)
            id_target_to_query = {}
            qid = 1
            tid = 1
            for i, (qaa, taa) in enumerate(zip(qaln_seq, taln_seq)):
                if taa != "-":
                    tid_ori = tid + tstart - 1
                    if tid_ori > t_key_sites_max:
                        break
                    if qaa != "-":
                        id_target_to_query[tid_ori] = qid + qstart - 1
                        qid += 1
                    tid += 1
                    continue
                if qaa != "-":
                    qid += 1
            curr_query_key_sites = [
                id_target_to_query.get(i, "unknown") for i in t_key_sites
            ]
            curr_qid_to_features = [
                f"{id_target_to_query[k]}: {v}"
                for k, v in curr_tid_to_features.items()
                if k in id_target_to_query
            ]
        query_key_sites.append(curr_query_key_sites)
        query_key_sites_features.append(" / ".join(curr_qid_to_features))
    search_result["query_key_sites"] = query_key_sites
    search_result["query_key_sites_features"] = query_key_sites_features
    search_result["query_key_sites"] = search_result["query_key_sites"].apply(
        lambda x: ";".join(map(str, x))
    )
    search_result.to_csv(
        str(Path(search_results_dir, "search_results.csv")), index=False
    )

    keysites = (
        search_result[["query", "query_key_sites"]]
        .groupby("query")
        .agg(agg_key_sites)
        .reset_index()
    )
    keysites["query_key_sites"] = keysites["query_key_sites"].apply(
        lambda x: ";".join(sorted(set([i for i in x.split(";") if i != "unknown"]), key=lambda x: int(x)))
    )
    keysites.to_csv(str(Path(search_results_dir, "key_sites.csv")), index=False)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Search key sites")
    parser.add_argument(
        "-i",
        "--docking_results_path",
        type=str,
        required=True,
        help="path to docking results csv file",
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=None, help="path to output dir"
    )
    parser.add_argument(
        "-n",
        "--n_workers",
        type=int,
        default=os.cpu_count(),
        help="number of workers",
    )
    parser.add_argument(
        "-e",
        "--multitasking_engine",
        type=str,
        default="process",
        choices=["process", "thread"],
        help="multitasking engine",
    )
    args = parser.parse_args()
    search_key_sites(**vars(args))


if __name__ == "__main__":
    main()
