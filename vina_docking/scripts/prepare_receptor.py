from pathlib import Path
import os
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool
import subprocess
from functools import partial
import warnings

warnings.filterwarnings("ignore")
BIN_LIB = str(Path(Path(__file__).parent.parent, "bin"))
PREPARE_RECEPTOR = str(Path(BIN_LIB, "ADFRsuit/ADFRsuite-1.1dev/bin/prepare_receptor"))


def run_prepare_receptor(pdb_filepath, output_dir):
    outfile = f"{output_dir}/{Path(pdb_filepath).stem}.pdbqt"
    if not os.path.exists(outfile):
        subprocess.run(
            [PREPARE_RECEPTOR, "-r", pdb_filepath, "-o", outfile, "-A", "hydrogens"],
            stdout=subprocess.PIPE,
        )


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-p", "--pdb_dir", type=str, default=None)
    parser.add_argument("-o", "--output_dir", type=str, default=None)
    args = parser.parse_args()
    pdb_files = glob(f"{args.pdb_dir}/*.pdb")
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    n_workers = os.cpu_count()
    pfunc = partial(run_prepare_receptor, output_dir=args.output_dir)
    with Pool(n_workers) as p:
        list(tqdm(p.imap(pfunc, pdb_files), total=len(pdb_files)))


if __name__ == "__main__":
    main()
