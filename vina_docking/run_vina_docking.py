#!/usr/bin/env python
import os
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool
from tempfile import TemporaryDirectory
import subprocess
import warnings
from pathlib import Path
import shutil
import pandas as pd
from functools import partial
from pprint import pprint
import sys
from scripts.prepare_receptor import ADFRSUITE
if ADFRSUITE not in sys.path:
    sys.path.append(ADFRSUITE)
from scripts import (
    run_docking,
    run_prepare_ligand,
    run_prepare_receptor,
    get_docking_details,
)
warnings.filterwarnings("ignore")


def to_sdf(ligand_path, out_path):
    """将ligand文件转换为3d sdf文件

    Args:
        ligand_path (str): 输入的ligand文件路径
        out_path (str): 输出的sdf文件路径, 需要以.sdf结尾
    """
    subprocess.run(
        ["obabel", ligand_path, "-O", out_path, "--gen3d"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def run_prepare_ligand_pipeline(ligand_path, output_dir):
    with TemporaryDirectory() as td:
        to_sdf(ligand_path=ligand_path, out_path=f"{td}/{Path(ligand_path).stem}.sdf")
        run_prepare_ligand(
            ligand_path=f"{td}/{Path(ligand_path).stem}.sdf", output_dir=output_dir
        )


def docking_pfunc(input_vals):
    """主要的docking函数, 用于多进程

    Args:
        input_vals (_type_): 包含receptor文件路径, ligand文件路径, docking_results_dir文件夹路径的list或者tuple

    Returns:
        tuple: 包含receptor名称, ligand名称, affinity, docking_result_dir文件夹路径的tuple
    """
    return run_docking(*input_vals)


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(add_help=True)
    parser.add_argument(
        "-l",
        "--ligand_dir",
        type=str,
        required=True,
        help="ligand文件夹路径, 支持mol, mol2, sdf格式",
    )
    parser.add_argument(
        "-r",
        "--receptor_dir",
        type=str,
        required=True,
        help="receptor文件夹路径, 支持cif(mmcif), pdb, pdbqt格式",
    )
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="输出文件夹路径")
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
    pprint(vars(args))
    ligand_files = glob(f"{args.ligand_dir}/*")
    receport_files = glob(f"{args.receptor_dir}/*")
    prepared_ligand_dir = f"{args.output_dir}/ligands_pdbqt"
    prepared_receptor_dir = f"{args.output_dir}/receptors_pdbqt"
    docking_results_dir = f"{args.output_dir}/docking_results"
    Path(prepared_ligand_dir).mkdir(parents=True, exist_ok=True)
    Path(prepared_receptor_dir).mkdir(parents=True, exist_ok=True)
    Path(docking_results_dir).mkdir(parents=True, exist_ok=True)
    with Pool(args.n_workers) as p:
        pfunc = partial(run_prepare_ligand_pipeline, output_dir=prepared_ligand_dir)
        list(
            tqdm(
                p.imap(pfunc, ligand_files),
                total=len(ligand_files),
                desc="prepare ligands",
            )
        )
    with Pool(args.n_workers) as p:
        pfunc = partial(run_prepare_receptor, output_dir=prepared_receptor_dir)
        list(
            tqdm(
                p.imap(pfunc, receport_files),
                total=len(receport_files),
                desc="prepare receptors",
            )
        )
    docking_inputs = [
        [receptor, ligand, docking_results_dir]
        for ligand in glob(f"{prepared_ligand_dir}/*")
        for receptor in glob(f"{prepared_receptor_dir}/*")
    ]
    with Pool(args.n_workers) as p:
        docking_results = list(
            tqdm(
                p.imap(docking_pfunc, docking_inputs),
                total=len(docking_inputs),
                desc="docking",
            )
        )
    result_df = pd.DataFrame(
        docking_results,
        columns=["receptor", "ligand", "affinity", "docking_result_dir"],
    )
    result_save_path = str(Path(args.output_dir, "docking_results.csv"))
    result_df.to_csv(result_save_path, index=None)
    get_docking_details(
        **{
            **vars(args),
            **dict(
                docking_results_path=result_save_path,
                out_filepath=str(Path(args.output_dir, "docking_results_details.csv")),
            ),
        }
    )


if __name__ == "__main__":
    main()
