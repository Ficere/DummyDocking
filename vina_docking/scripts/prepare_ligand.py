import os
from pathlib import Path
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
import subprocess
import warnings

warnings.filterwarnings("ignore")


def run_prepare_ligand(ligand_path, output_dir):
    out_file = f"{output_dir}/{Path(ligand_path).stem}.pdbqt"
    subprocess.run(
        [
            "mk_prepare_ligand.py",
            "-i",
            ligand_path,
            "-o",
            out_file,
            "--merge_these_atom_types",
        ],
        stdout=subprocess.PIPE,
    )


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-l", "--ligand_dir", type=str, default=None)
    parser.add_argument("-o", "--output_dir", type=str, default=None)
    args = parser.parse_args()
    ligand_files = glob(f"{args.ligand_dir}/*")
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    n_workers = os.cpu_count()
    with Pool(n_workers) as p:
        pfunc = partial(run_prepare_ligand, output_dir=args.output_dir)
        list(tqdm(p.imap(pfunc, ligand_files), total=len(ligand_files)))


if __name__ == "__main__":
    main()
