import subprocess
from pathlib import Path
import os

#python /home/guolj/ws/github/DummyDocking/vina_docking/run_vina_docking.py -r /home/guolj/ws/github/DummyDocking/example/receptor -l /home/guolj/ws/github/DummyDocking/example/ligand -o /home/guolj/ws/github/DummyDocking/example/output -k
def main():
    curr_dir = Path(__file__).parent
    subprocess.run(
        [
            "python",
            str(curr_dir / "vina_docking" / "run_vina_docking.py"),
            "-r",
            str(curr_dir / "example" / "receptor"),
            "-l",
            str(curr_dir / "example" / "ligand"),
            "-o",
            str(curr_dir / "example" / "output"),
            "-k",
        ],
    )

if __name__ == "__main__":
    main()