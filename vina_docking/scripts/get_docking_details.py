import numpy as np
from scipy.spatial import distance_matrix
import pandas as pd
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser
from pathlib import Path
from glob import glob
from tqdm import tqdm
import subprocess
from typing import Optional
import warnings
import os
from itertools import chain
import re
if __name__ == "__main__":
    from search_key_sites import search_key_sites
else:
    from .search_key_sites import search_key_sites

warnings.filterwarnings("ignore")


def get_docking_details(
    docking_results_path: str,
    out_filepath: str,
    n_workers: int = os.cpu_count(),
    multitasking_engine: str = "process",
    distance_threshold: float = 4.0,
    do_search_key_sites: bool = True,
    buffer: int = 5,
    search_key_sites_outdir: Optional[str] = None,
    **kwargs,
):
    docking_results = pd.read_csv(docking_results_path)
    out = []
    for receptor_name, ligand_name, curr_result_dir in tqdm(
        docking_results[
            ["receptor", "ligand", "docking_result_dir"]
        ].values,
        total=len(docking_results),
        desc="get docking details",
    ):
        curr_receptor = str(Path(curr_result_dir, "receptor.pdbqt"))
        curr_ligand = str(Path(curr_result_dir, "docking_result.pdbqt"))
        curr_log = str(Path(curr_result_dir, "vina_dock.log"))

        with open(curr_log, "r") as f:
            log = f.readlines()
            for i, v in enumerate(log):
                if v.startswith("mode"):
                    scores = list(
                        map(lambda x: re.split(" +", x.strip()), log[i + 3 :])
                    )
                    break
            else:
                raise ValueError("No mode found in log file")

        scores = pd.DataFrame(
            scores, columns=["mode", "affinity", "rmsd_lb", "rmsd_ub"]
        )

        parser = PDBParser()
        receptor_pdb = parser.get_structure("", curr_receptor)[0]
        ligand_pdb = parser.get_structure("", curr_ligand)

        receptor_resnames = np.array(
            [
                "".join(
                    [
                        seq1(atom.parent.resname),
                        "".join(map(str, atom.parent.id)).strip(),
                    ]
                )
                for atom in receptor_pdb.get_atoms()
            ]
        )
        receptor_resnames_to_id = {
            "".join(
                [
                    seq1(res.resname),
                    "".join(map(str, res.id)).strip(),
                ]
            ): i + 1
            for i, res in enumerate(receptor_pdb.get_residues())
        }
        receptor_coords = np.array(
            [atom.get_coord() for atom in receptor_pdb.get_atoms()]
        )

        for curr_ligand_model, affinity in zip(
            ligand_pdb, scores["affinity"].values
        ):
            curr_ligand_coords = np.array(
                [atom.get_coord() for atom in curr_ligand_model.get_atoms()]
            )
            dis_mat = distance_matrix(receptor_coords, curr_ligand_coords)
            interact_residues = np.unique(
                receptor_resnames[dis_mat.min(axis=1) < distance_threshold]
            )
            out.append(
                [
                    receptor_name,
                    ligand_name,
                    curr_result_dir,
                    ";".join(interact_residues),
                    ";".join(
                        [str(receptor_resnames_to_id[i]) for i in interact_residues]
                    ),
                    affinity,
                ]
            )
    out = pd.DataFrame(
        out,
        columns=[
            "receptor",
            "ligand",
            "docking_result_dir",
            "interact_residues",
            "interact_residue_ids",
            "affinity",
        ],
    )
    if do_search_key_sites:
        print("Searching key sites...")
        search_key_sites_outdir = (
            str(Path(Path(docking_results_path).parent, "search_results"))
            if search_key_sites_outdir is None
            else str(search_key_sites_outdir)
        )
        search_key_sites(docking_results_path, search_key_sites_outdir, n_workers, multitasking_engine)
        print("Done searching key sites")
        print("Merging key sites with docking results...")
        key_sites = pd.read_csv(str(Path(search_key_sites_outdir, "key_sites.csv")))
        query_to_key_sites = {
            k: list(map(int, v.split(";")))
            for k, v in key_sites[["query", "query_key_sites"]].values
        }
        key_sites_involved = []
        for receptor_name, interact_residue_ids in out[["receptor", "interact_residue_ids"]].values:
            interact_residue_ids_with_buffer = [list(range(i - buffer, i + buffer + 1)) for i in map(int, interact_residue_ids.split(";"))]
            interact_residue_ids_with_buffer = set(chain(*interact_residue_ids_with_buffer))
            key_sites_involved.append(any([i in query_to_key_sites.get(receptor_name, []) for i in interact_residue_ids_with_buffer]))
        out["key_sites_involved"] = key_sites_involved
    out.to_csv(out_filepath, index=False)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--docking_results_path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out_filepath",
        type=str,
        required=True,
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
    parser.add_argument(
        "-d",
        "--distance_threshold",
        type=float,
        default=4.0,
        help="Distance threshold (in terms of angstrom) for identifying interacting residues",
    )
    parser.add_argument(
        "-k",
        "--do_search_key_sites",
        action="store_true",
        help="Search key sites",
    )
    parser.add_argument(
        "-b",
        "--buffer",
        type=int,
        default=5,
        help="Buffer for identifying key sites. Residues within buffer distance (in terms of residue id) of interacting residues will be considered as involved in the interaction",
    )
    parser.add_argument(
        "-ko",
        "--search_key_sites_outdir",
        type=str,
        default=None,
        help="Path to output dir for search key sites",
    )
    args = parser.parse_args()
    get_docking_details(**vars(args))
    


if __name__ == "__main__":
    main()
