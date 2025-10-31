import shutil
from pathlib import Path
import subprocess 
import math
import os
from typing import Tuple, Union, List, Dict, Any, Optional


def which_or_raise(bin_name: str):
    if shutil.which(bin_name) is None:
        raise RuntimeError(f"Required tool not found on PATH: {bin_name}")
    
def run(cmd: Union[List[str], str], **kw) -> None:
    print("[sh]", cmd if isinstance(cmd, str) else " ".join(cmd))
    subprocess.run(cmd, check=True, shell=isinstance(cmd, str), **kw)

def ensure_fasta_faidx(ref_fa: Path) -> None:
    which_or_raise("samtools")
    if not ref_fa.exists():
        raise FileNotFoundError(f"FASTA not found: {ref_fa}")
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if not fai.exists():
        run(["samtools", "faidx", str(ref_fa)])

def ensure_vcf_bgzip_tabix(vcf_path: Path) -> Path:
    which_or_raise("bgzip"); which_or_raise("tabix")
    if vcf_path.suffixes[-2:] == [".vcf", ".gz"]:
        vcfgz = vcf_path
    elif vcf_path.suffix == ".vcf":
        vcfgz = vcf_path.with_suffix(vcf_path.suffix + ".gz")
        if not vcfgz.exists():
            run(f"bgzip -c {vcf_path} > {vcfgz}")
    elif vcf_path.suffix == ".gz":
        vcfgz = vcf_path
    else:
        raise ValueError(f"Unsupported VCF path: {vcf_path}")
    tbi = Path(str(vcfgz) + ".tbi")
    if not tbi.exists():
        run(["tabix", "-f", "-p", "vcf", str(vcfgz)])
    return vcfgz

def prepare_inputs(ref_fa: Path, vcf: Path) -> Tuple[Path, Path]:
    ensure_fasta_faidx(ref_fa)
    vcfgz = ensure_vcf_bgzip_tabix(vcf)
    return ref_fa, vcfgz

def _need(bin_name: str):
    if shutil.which(bin_name) is None:
        raise RuntimeError(f"Tool not found on PATH: {bin_name}")

def _pairs_for(depth: float, L: int, read_len: int) -> int:
    return max(1, math.ceil(depth * L / (2.0 * read_len)))

def _mkfifo(p: Path):
    if p.exists():
        try: os.remove(p)
        except FileNotFoundError: pass
    os.mkfifo(p)
    return p

def _spawn(cmd: str) -> subprocess.Popen:
    """Spawn a shell command (string) and return Popen."""
    return subprocess.Popen(cmd, shell=True, executable="/bin/bash")
import shlex
import time

def fix_ref_by_vcf_ref(vcfgz: Path, ref_template: Path, ref_fixed: Path) -> Path:
    """
    用 VCF 的 REF (-H R) 覆盖模板参考（可为全A），生成与 VCF 完全一致的 ref_fixed.fa。
    只为 ref_fixed 建索引（.fai）。
    """
    which_or_raise("bcftools"); which_or_raise("samtools")

    ref_fixed.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(ref_fixed, "w") as fout:
            subprocess.run(
                ["bcftools", "consensus", "-H", "R", "-f", str(ref_template), str(vcfgz)],
                check=True, stdout=fout, stderr=subprocess.PIPE, text=True
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"[bcftools consensus -H R] failed.\n"
            f"CMD: bcftools consensus -H R -f {ref_template} {vcfgz}\n"
            f"STDERR:\n{e.stderr}"
        ) from e

    # 只给修正后的参考建索引
    subprocess.run(["samtools", "faidx", str(ref_fixed)], check=True)
    return ref_fixed


def simulate_reads_streaming(
    sample: str,
    *,
    ref_fa: Path,          # reference FASTA (indexed with samtools faidx)
    vcf_gz: Path,          # bgzip-compressed VCF with tabix index
    out_dir: Path,         # output directory for {sample}_R1.fq.gz / {sample}_R2.fq.gz
    L: int,                # target region length (used to compute -N)
    depth: float = 10.0,   # target coverage
    read_len: int = 150,   # read length for R1/R2
    err: float = 0.001,    # base error rate
    insert_mean: int = 350,# fragment length mean
    insert_sd: int = 35    # fragment length stddev
) -> Tuple[Path, Path]:
    """
    Generate paired-end reads with wgsim without writing hapFASTA or raw FASTQ to disk.
    Pipeline per hap:
      bcftools consensus (-H 1/2) -> hap_fifo -> wgsim -> (r1_fifo,r2_fifo) -> gzip -> {sample}_R1/R2.fq.gz
    hap1 uses truncate write (>), hap2 appends (>>) into the same gz files.
    """
    # dependencies
    for t in ("bcftools", "wgsim", "gzip"):
        _need(t)

    # outputs
    out_dir.mkdir(parents=True, exist_ok=True)
    out_r1 = out_dir / f"{sample}_R1.fq.gz"
    out_r2 = out_dir / f"{sample}_R2.fq.gz"
    if out_r1.exists():
        out_r1.unlink()
    if out_r2.exists():
        out_r2.unlink()

    ref_fixed = out_dir / ".ref.fixed.fa"
    fix_ref_by_vcf_ref(vcf_gz, ref_fa, ref_fixed)

    # read-pair plan across two haplotypes
    total_pairs = _pairs_for(depth, L, read_len)
    half = max(1, total_pairs // 2)
    plan = [(1, half, False), (2, total_pairs - half, True)]  # (hap, npairs, append?)

    for hap, npairs, append in plan:
        # FIFOs
        hap_fifo = _mkfifo(out_dir / f".{sample}.hap{hap}.fa.fifo")
        r1_fifo  = _mkfifo(out_dir / f".{sample}.hap{hap}.R1.fq.fifo")
        r2_fifo  = _mkfifo(out_dir / f".{sample}.hap{hap}.R2.fq.fifo")

        try:
            # gzip sinks first (listen on R1/R2 FIFOs)
            redir = ">>" if append else ">"
            gz1 = _spawn(f"gzip -c < {shlex.quote(str(r1_fifo))} {redir} {shlex.quote(str(out_r1))}")
            gz2 = _spawn(f"gzip -c < {shlex.quote(str(r2_fifo))} {redir} {shlex.quote(str(out_r2))}")
            time.sleep(0.05)

            # bcftools consensus -> hap_fifo (no hapFASTA on disk)
            cons = _spawn(
                " ".join([
                    "bcftools consensus",
                    f"-f {shlex.quote(str(ref_fixed))}",
                    f"-s {shlex.quote(sample)}",
                    f"-H {hap}",
                    shlex.quote(str(vcf_gz)),
                    f"> {shlex.quote(str(hap_fifo))}",
                ])
            )
            time.sleep(0.05)

            # wgsim: hap_fifo -> r1_fifo/r2_fifo (raw FASTQ to gzip streams)
            wgs = _spawn(
                " ".join([
                    "wgsim",
                    f"-e {err}",
                    f"-d {insert_mean}",
                    f"-s {insert_sd}",
                    f"-N {npairs}",
                    f"-1 {read_len}",
                    f"-2 {read_len}",
                    shlex.quote(str(hap_fifo)),
                    shlex.quote(str(r1_fifo)),
                    shlex.quote(str(r2_fifo)),
                ])
            )

            # wait for this hap's pipeline
            wgs.wait()
            gz1.wait(); gz2.wait()
            cons.wait()

        finally:
            # clean FIFOs
            for p in (hap_fifo, r1_fifo, r2_fifo):
                try:
                    os.remove(p)
                except FileNotFoundError:
                    pass

    return out_r1, out_r2