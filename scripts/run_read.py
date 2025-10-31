# scripts/run_read.py
import gzip
import yaml
from pathlib import Path
# 顶部增加
from typing import List


from src.read import prepare_inputs   # 负责：ref.fa.fai、VCF→.vcf.gz+.tbi、contig检查
from src.read import simulate_reads_streaming


def _vcf_sample_ids(vcf_path: Path) -> List[str]:
    """
    读取 VCF 头，返回样本名列表（#CHROM 行第9列起）。
    支持 .vcf 和 .vcf.gz
    """
    opener = gzip.open if vcf_path.suffix.endswith("gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith("#CHROM"):
                cols = ln.strip().split("\t")
                return cols[9:]  # 9列起是样本
    return []


def run_read():
    cfg = yaml.safe_load(Path("config/default.yaml").read_text())

    outdir = Path(cfg["project"]["outdir"])
    sim_dir = outdir / "simulate"     # 和你的 run_simulate 输出保持一致
    reads_dir = outdir / "reads"      # 本步输出目录

    # 输入路径（与 run_simulate 的输出约定一致）
    ref_fa = sim_dir / "ref.fa"
    vcf_in = sim_dir / "truth.vcf"    # 可能是 .vcf 或 .vcf.gz，prepare_inputs 会处理

    # 预处理：确保 ref.fa.fai、把 VCF 变成 .vcf.gz 并建 .tbi
    ref_fa, vcf_gz = prepare_inputs(ref_fa, vcf_in)

    # 读取 wgsim 参数（带默认值，缺了也能跑）
    msprime_l = int(cfg["msprime"]["l"])
    reads_cfg  = cfg.get("reads", {})
    depth      = float(reads_cfg.get("depth", 10.0))
    read_len   = int(reads_cfg.get("read_len", 150))
    err        = float(reads_cfg.get("err", 0.001))
    insert_mu  = int(reads_cfg.get("insert_mean", 350))
    insert_sd  = int(reads_cfg.get("insert_sd", 35))
    # L 优先用 reads.L，否则回退到 msprime.l
    L          = int(reads_cfg.get("L", msprime_l))

    # 样本列表：优先用配置里的 reads.samples；否则从 VCF 头部读
    samples = reads_cfg.get("samples")
    if samples is None:
        samples = _vcf_sample_ids(vcf_gz)
    max_samples = reads_cfg.get("max_samples")
    if max_samples is not None:
        samples = samples[: int(max_samples)]

    if not samples:
        raise RuntimeError("No samples found in VCF and no 'reads.samples' provided in config.")

    reads_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Generating reads for {len(samples)} samples → {reads_dir}")
    for i, sample in enumerate(samples, 1):
        print(f"[{i}/{len(samples)}] {sample} …", end="", flush=True)
        r1, r2 = simulate_reads_streaming(
            sample=sample,
            ref_fa=ref_fa,
            vcf_gz=vcf_gz,
            out_dir=reads_dir,
            L=L,
            depth=depth,
            read_len=read_len,
            err=err,
            insert_mean=insert_mu,
            insert_sd=insert_sd,
        )
        print(f" done ({r1.name}, {r2.name})")

    print("[OK] Read simulation complete.")


if __name__ == "__main__":
    run_read()
