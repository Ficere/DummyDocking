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
conda activate vina_docking
conda install -yc conda-forge mamba
mamba install -yc conda-forge rdkit
mamba install -yc bioconda mmseqs2
pip install panda meeko openpyxl scipy tqdm biopython biopandas
```

# 使用
```
usage: run_vina_docking.py [-h] -l LIGAND_DIR -r RECEPTOR_DIR -o OUTPUT_DIR [-n N_WORKERS] [-e {process,thread}] [-d DISTANCE_THRESHOLD] [-k] [-b BUFFER] [-ko SEARCH_KEY_SITES_OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -l LIGAND_DIR, --ligand_dir LIGAND_DIR
                        ligand文件夹路径, 支持mol, mol2, sdf格式
  -r RECEPTOR_DIR, --receptor_dir RECEPTOR_DIR
                        receptor文件夹路径, 支持cif(mmcif), pdb, pdbqt格式
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        输出文件夹路径
  -n N_WORKERS, --n_workers N_WORKERS
                        工作进程数量
  -e {process,thread}, --multitasking_engine {process,thread}
                        multitasking所用的引擎
  -d DISTANCE_THRESHOLD, --distance_threshold DISTANCE_THRESHOLD
                        用于识别相互作用残基的距离阈值（以埃为单位）
  -k, --do_search_key_sites
                        在Swiss-Prot上搜索关键位点
  -b BUFFER, --buffer BUFFER
                        用于识别关键位点的缓冲区。与小分子相互作用残基的距离在缓冲区距离（以残基序号为单位）内的残基将被视为参与相互作用
  -ko SEARCH_KEY_SITES_OUTDIR, --search_key_sites_outdir SEARCH_KEY_SITES_OUTDIR
                        搜索关键位点的输出目录的路径
```

我在example中为你准备了一些测试数据, 你可以在example/output中看到示例输出

当然，你也可以这样进行测试

```bash
# 请确保vina_docking/run_vina_docking.py是可执行文件
# 你可以简单地chmod 755 vina_docking/run_vina_docking.py实现
vina_docking/run_vina_docking.py -l example/ligand -r example/receptor -o example/output -k
# 或者你也可以自己指定python解释器位置
python vina_docking/run_vina_docking.py -l example/ligand -r example/receptor -o example/output -k
```