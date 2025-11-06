# scripts/run_read.py
import gzip
import yaml
from pathlib import Path
# 顶部增加
from typing import List

from src.read import SequencingRead  # 负责：从 VCF 构建 haplotypes、写 FASTA


def run_read():
    cfg = yaml.safe_load(Path("config/default.yaml").read_text())

    SequencingRead(
        vcf_path=str(cfg["project"]["simulate_vcfs"]),
        L=cfg["msprime"]["l"],
        fasta_path=str(cfg["project"]["read_all"]),
    ).write_all_sample_haplotypes_fasta()

    print("[OK] Read simulation complete.")


if __name__ == "__main__":
    run_read()
