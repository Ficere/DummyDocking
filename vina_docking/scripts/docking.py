import os
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool
import subprocess
import warnings
from rdkit import Chem
from pathlib import Path
import shutil
import re

warnings.filterwarnings("ignore")

BIN_LIB = str(Path(Path(__file__).parent.parent, "bin"))
AUTOGRID = str(Path(BIN_LIB, "ADFRsuit/ADFRsuite-1.1dev/bin/autogrid4"))
PYTHON2SH = str(Path(BIN_LIB, "ADFRsuit/ADFRsuite-1.1dev/bin/pythonsh"))
PREPARE_GDF = str(
    Path(BIN_LIB, "AutoDock-Vina/example/autodock_scripts/prepare_gpf.py")
)
MAPWATER = str(Path(BIN_LIB, "AutoDock-Vina/example/autodock_scripts/mapwater.py"))
VINA = str(Path(BIN_LIB, "vina_1.2.3_linux_x86_64"))
DRYPY = str(Path(BIN_LIB, "AutoDock-Vina/example/autodock_scripts/dry.py"))

def best_affinity(log_file):
    with open(log_file, "r") as f:
        res = [
            float(re.split(" +", line.strip())[1])
            for line in f
            if re.search("^ +\d", line)
        ]
        cur_low = min(res) if res else None
    return cur_low

def run_docking(receptor_path, ligand_path, output_dir):
    workdir = Path(output_dir, f"{Path(receptor_path).stem};;{Path(ligand_path).stem}")
    Path(workdir).mkdir(parents=True, exist_ok=True)
    shutil.copy(receptor_path, Path(workdir, "receptor.pdbqt"))
    shutil.copy(ligand_path, Path(workdir, "ligand.pdbqt"))
    subprocess.run(
        [
            PYTHON2SH,
            PREPARE_GDF,
            "-l",
            "ligand.pdbqt",
            "-r",
            "receptor.pdbqt",
            "-o",
            "gpf_result.gpf",
        ],
        cwd=workdir,
        stdout=open(Path(workdir, "prepare_gdf.log"), "w"),
        stderr=open(Path(workdir, "prepare_gdf.err"), "w"),
    )
    subprocess.run(
        [
            AUTOGRID,
            "-p",
            "gpf_result.gpf",
            "-l",
            "autogrid_result.glg",
        ],
        cwd=workdir,
        stdout=open(Path(workdir, "autogrid.log"), "w"),
        stderr=open(Path(workdir, "autogrid.err"), "w"),
    )
    subprocess.run(
        [
            VINA,
            "--ligand",
            "ligand.pdbqt",
            "--maps",
            "receptor",
            "--scoring",
            "ad4",
            "--exhaustiveness",
            "8",
            "--out",
            "docking_result.pdbqt",
        ],
        cwd=workdir,
        stdout=open(Path(workdir, "vina_dock.log"), "w"),
        stderr=open(Path(workdir, "vina_dock.err"), "w"),
    )
    affinity = best_affinity(Path(workdir, "vina_dock.log"))
    return Path(receptor_path).stem, Path(ligand_path).stem, affinity, workdir

def docking_pfunc(input_vals):
    run_docking(*input_vals)

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-r", "--receptor_pdbqt_dir", type=str, default=None)
    parser.add_argument("-l", "--ligand_pdbqt_dir", type=str, default=None)
    parser.add_argument("-o", "--output_dir", type=str, default=None)
    parser.add_argument("-n", "--n_cpus", type=int, default=os.cpu_count())
    args = parser.parse_args()
    ligand_files = glob(f"{args.ligand_pdbqt_dir}/*")
    receptor_files = glob(f"{args.receptor_pdbqt_dir}/*")
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    pairs = [(i, j, args.output_dir) for i in receptor_files for j in ligand_files]
    with Pool(args.n_cpus) as p:
        list(tqdm(p.imap(docking_pfunc, pairs), total=len(pairs)))


if __name__ == "__main__":
    main()
