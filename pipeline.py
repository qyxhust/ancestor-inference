# pipeline.py —— 极简主控：调用你已有的 run_simulate() 跑 msprime
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
# 让 Python 能找到 scripts/ 与 src/（src 采用包式导入 from src.xxx import ...）
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

from scripts.run_simulate import run_simulate  # 你的函数，内部已内置 cfg_path
from scripts.run_read import run_read          # 你的函数，内部已内置 cfg_path


def main():
    run_simulate()
    print("[OK] Pipeline complete.")
    run_read()
    print("[OK] Read simulation complete.")

if __name__ == "__main__":
    main()
