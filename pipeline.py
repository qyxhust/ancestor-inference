# pipeline.py —— 极简主控：调用你已有的 run_simulate() 跑 msprime
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
# 让 Python 能找到 scripts/ 与 src/（src 采用包式导入 from src.xxx import ...）
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

from scripts.run_simulate import run_simulate  # 你的函数，内部已内置 cfg_path
from scripts.run_read import run_read          # 你的函数，内部已内置 cfg_path
from src.read import stratified_split_samples_pd  # 负责：从 VCF 构建 haplotypes、写 FASTA


def main():
    # 1) 运行模拟，生成 VCF + 标签
    run_simulate()
    print("[OK] Pipeline complete.")

    # 2) 分层划分样本
    stratified_split_samples_pd(
        meta_path="data/simulate/labels.tsv",
        ref_out="data/ref_samples.csv",
        test_out="data/test_samples.csv",
        ref_ratio=0.75,
        seed=42,
        sep="\t",
    )

    run_read()
    print("[OK] Read simulation complete.")

if __name__ == "__main__":
    main()
