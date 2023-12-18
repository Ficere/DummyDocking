![License](https://img.shields.io/badge/license-MIT-yellowgreen)  ![Language](https://img.shields.io/badge/language-python-blue)

<a href="./readme.chs.md">简体中文</a> | <a href="./readme.md">English</a>

# Basic Vina Docking
给我这样的笨比使用的分子对接流程, 包含了:

1. 自动配体和受体准备

2. 自动对接并获取分数（目前只支持basic docking）, 然后输出到一个csv文件中

# 环境安装
```bash
# 来个新的环境
conda create -yn vina_docking python=3.10
conda install -yc conda-forge mamba
mamba install -yc conda-forge rdkit openbabel
pip install pandas meeko openpyxl scipy
```

# 使用
```
usage: vina_docking/run_vina_docking.py [-h] -l LIGAND_DIR -r RECEPTOR_DIR -o OUTPUT_DIR [-n N_WORKERS]

optional arguments:
  -h, --help            show this help message and exit
  -l LIGAND_DIR, --ligand_dir LIGAND_DIR
                        ligand文件夹路径, 支持mol, mol2, sdf格式
  -r RECEPTOR_DIR, --receptor_dir RECEPTOR_DIR
                        receptor文件夹路径, 支持cif(mmcif), pdb, pdbqt格式
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        输出文件夹路径
  -n N_WORKERS, --n_workers N_WORKERS
```

我在example中为你准备了一些测试数据, 你可以在example/output中看到示例输出

当然，你也可以这样进行测试

```bash
python vina_docking/run_vina_docking.py -l example/ligand -r example/receptor -o example/output
```
