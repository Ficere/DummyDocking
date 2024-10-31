![License](https://img.shields.io/badge/license-MIT-yellowgreen)  ![Language](https://img.shields.io/badge/language-python-blue)

<a href="./readme.chs.md">简体中文</a> | <a href="./readme.md">English</a>

# Basic Vina Docking
Here is a molecular docking procedure designed for dummies like me, including:

1. Automatic preparation of ligands and receptors.

2. Automated docking with score retrieval (only basic docking implemented yet). outputing to a CSV file.

3. Retrieve interacting residues based on docking results and output them to another CSV file.

4. (Optional) Search the SProt database based on the sequence to obtain key residues (known binding sites, etc.), and determine whether the interaction pocket protects these key residues.

# Installing Dependencies
```bash
# Set up a new environment
conda create -yn vina_docking python=3.10
conda activate vina_docking
conda install -yc conda-forge mamba
mamba install -yc conda-forge openbabel
mamba install -yc bioconda mmseqs2
pip install pandas meeko openpyxl scipy tqdm biopython multitasking retry lxml rdkit
git clone https://github.com/Ficere/DummyDocking.git 
cd DummyDocking/vina_docking/bin/ADFRsuit 
chmod 777 install.sh
./install.sh
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
                        number of workers
  -e {process,thread}, --multitasking_engine {process,thread}
                        multitasking engine
  -d DISTANCE_THRESHOLD, --distance_threshold DISTANCE_THRESHOLD
                        Distance threshold (in terms of angstrom) for identifying interacting residues
  -k, --do_search_key_sites
                        Search key sites in Swiss-Prot
  -b BUFFER, --buffer BUFFER
                        Buffer for identifying key sites. Residues within buffer distance (in terms of residue id) of interacting residues will be considered as involved in the interaction
  -ko SEARCH_KEY_SITES_OUTDIR, --search_key_sites_outdir SEARCH_KEY_SITES_OUTDIR
                        Path to output dir for search key sites
```

I've prepared some test data for you in the example folder, and you can find sample outputs in the example/output directory.

Of course, you can also perform the test like this:

```bash
# Make sure vina_docking/run_vina_docking.py is executable
# You can easily achieve this by running `chmod 755 vina_docking/run_vina_docking.py`
# Then, it should be
vina_docking/run_vina_docking.py -l example/ligand -r example/receptor -o example/output
# Alternatively, you can specify the location of the Python interpreter yourself.
python vina_docking/run_vina_docking.py -l example/ligand -r example/receptor -o example/output -k
```
