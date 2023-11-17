![License](https://img.shields.io/badge/license-MIT-yellowgreen)  ![Language](https://img.shields.io/badge/language-python-blue)

<a href="./readme.chs.md">简体中文</a> | <a href="./readme.md">English</a>

# Basic Vina Docking
Here is a molecular docking procedure designed for dummies like me, including:

1. Automatic preparation of ligands and receptors.

2. Automated docking with score retrieval (only basic docking implemented yet). outputing to a CSV file.

# Installation 
```bash
# Set up a new environment
conda create -yn vina_docking python=3.10
conda install -yc conda-forge mamba
mamba install -yc conda-forge rdkit
pip install pandas meeko openpyxl scipy
```

# Usage
```
usage: vina_docking/run_vina_docking.py [-h] -l LIGAND_DIR -r RECEPTOR_DIR -o OUTPUT_DIR [-n N_WORKERS]

optional arguments:
  -h, --help            show this help message and exit
  -l LIGAND_DIR, --ligand_dir LIGAND_DIR
                        path to the ligand folder, supporting mol, mol2, sdf formats
  -r RECEPTOR_DIR, --receptor_dir RECEPTOR_DIR
                        path to the receptor folder, supporting cif(mmcif), pdb, pdbqt formats
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        path to the output folder
  -n N_WORKERS, --n_workers N_WORKERS
```

I've prepared some test data for you in the example folder, and you can find sample outputs in the example/output directory.

Of course, you can also perform the test like this:

```bash
python vina_docking/run_vina_docking.py -l example/ligand -r example/receptor -o example/output
```
