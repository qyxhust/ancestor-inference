import math
import subprocess
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
import textwrap


class SequencingRead:
    """
    一条龙功能：
      1. 生成一条全局随机 ATGC 序列（base_random）
      2. 用 VCF 的 REF 覆盖所有 SNP -> 全局 ref_array
      3. 对每个样本，在 ref_array 基础上按 GT 改 ALT -> hap1/hap2
      4. 75% reference 样本的 hap 合并写到一个总 FASTA（EM 参考库）
      5. 25% test 样本：每个样本一个目录，写 sample.fa + wgsim 双端 reads
    """

    def __init__(
        self,
        vcf_path: str,
        L: int,
        read_root: str,
        chrom: str = "1",
        test_table: Optional[str] = None,
        sample_col: str = "sample",
        pop_col: str = "population",
        sep: str = ",",
        seed: int = 42,
    ):
        """
        vcf_path : msprime 生成的 VCF
        L        : 序列长度（>= VCF 中最大坐标）
        read_root: 输出的根目录（下面会有 reads/<sample>/ 等）
        chrom    : 处理的染色体
        test_table: 测试样本列表（TSV/CSV, 至少含 sample_col + pop_col）
                    -> 只对这些样本做测序（25%）
        seed     : 随机种子，用于生成基础随机序列
        """
        self.vcf_path = Path(vcf_path)
        self.L = int(L)
        self.read_root = Path(read_root)  # 一般用 "data/reads"
        self.chrom = chrom

        self.vcf = pysam.VariantFile(str(self.vcf_path))
        self.all_samples: List[str] = list(self.vcf.header.samples)

        # 载入 test/reference 划分
        if test_table is not None:
            meta = pd.read_csv(test_table, sep=sep)
            meta = meta[[sample_col, pop_col]].drop_duplicates()
            meta = meta.rename(columns={sample_col: "sample", pop_col: "population"})
            self.meta = meta
            self.test_samples = set(meta["sample"])
        else:
            self.meta = pd.DataFrame({"sample": self.all_samples, "population": "UNKNOWN"})
            self.test_samples = set(self.all_samples)

        self.reference_samples = [s for s in self.all_samples if s not in self.test_samples]

        print(f"[INFO] VCF samples total = {len(self.all_samples)}")
        print(f"[INFO] reference samples = {len(self.reference_samples)}, "
              f"test samples = {len(self.test_samples)}")

        # ===== 核心三步：随机基础序列 -> REF 覆盖 =====
        self.rng = np.random.default_rng(seed)
        self.base_random_array = self._build_base_random_array()
        self.ref_array = self._build_ref_array_from_vcf()

    # =====================================================
    #  功能 1：生成一条随机 ATGC 序列（只生成一次）
    # =====================================================

    def _build_base_random_array(self) -> np.ndarray:
        """
        生成长度为 L 的随机 ATGC 序列：
        - 只在 __init__ 里调用一次
        - 后续所有 ref/hap 都基于同一条序列
        """
        bases = np.array(list("ATGC"))
        idx = self.rng.integers(0, 4, size=self.L)
        base_random = bases[idx]
        print("[INFO] base_random_array generated.")
        return base_random

    # =====================================================
    #  功能 2：用 VCF 中的 REF 覆盖 SNP -> 全局 ref_array
    # =====================================================

    def _build_ref_array_from_vcf(self) -> np.ndarray:
        """
        在 base_random_array 基础上，用 VCF 中的 REF 覆盖所有 SNP 位点，
        得到全局参考序列 ref_array。
        """
        ref_array = self.base_random_array.copy()
        vcf = pysam.VariantFile(str(self.vcf_path))

        try:
            iterator = vcf.fetch(self.chrom)
        except (ValueError, TypeError):
            iterator = vcf

        for rec in iterator:
            if rec.contig != self.chrom:
                continue
            if len(rec.ref) != 1:
                continue

            pos = rec.pos - 1  # 1-based -> 0-based
            if 0 <= pos < self.L:
                ref_array[pos] = rec.ref

        print("[INFO] ref_array built from VCF REF.")
        return ref_array

    # =====================================================
    #  功能 3：对单个样本生成 hap1 / hap2（基于 ref_array）
    # =====================================================

    def _build_sample_haplotypes(self, sample: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        对单个样本：
          - 从 ref_array 拷贝两份
          - 按该样本的 GT，把 SNP 位点替换成 ALT
        返回：hap1_array, hap2_array
        """
        if sample not in self.all_samples:
            raise ValueError(f"Sample {sample} not found in VCF.")

        hap1 = self.ref_array.copy()
        hap2 = self.ref_array.copy()

        vcf = pysam.VariantFile(str(self.vcf_path))
        try:
            iterator = vcf.fetch(self.chrom)
        except (ValueError, TypeError):
            iterator = vcf

        for rec in iterator:
            if rec.contig != self.chrom:
                continue
            if len(rec.ref) != 1:
                continue

            pos = rec.pos - 1
            if not (0 <= pos < self.L):
                continue

            alts = rec.alts
            if not alts:
                continue

            sample_data = rec.samples[sample]
            gt = sample_data.get("GT", None)  # (0,1),(1,1)...
            if gt is None:
                continue

            for hap_idx, allele_index in enumerate(gt):
                if allele_index is None or allele_index == 0:
                    continue  # REF
                if allele_index - 1 >= len(alts):
                    continue

                alt_base = alts[allele_index - 1]
                if len(alt_base) != 1:
                    continue  # 跳过 indel

                if hap_idx == 0:
                    hap1[pos] = alt_base
                else:
                    hap2[pos] = alt_base

        return hap1, hap2

    # =====================================================
    #  写 FASTA 的小工具
    # =====================================================

    @staticmethod
    def _write_two_haps_to_fasta(
        sample: str,
        hap1: np.ndarray,
        hap2: np.ndarray,
        out_fa: Path,
        line_width: int = 100,
    ):
        seq1 = "".join(hap1)
        seq2 = "".join(hap2)
        out_fa.parent.mkdir(parents=True, exist_ok=True)

        with open(out_fa, "w") as f:
            f.write(f">{sample}_hap1\n")
            for chunk in textwrap.wrap(seq1, line_width):
                f.write(chunk + "\n")

            f.write(f">{sample}_hap2\n")
            for chunk in textwrap.wrap(seq2, line_width):
                f.write(chunk + "\n")

    # =====================================================
    #  4a. 写 75% reference：合并到一个总 FASTA（EM 参考库）
    # =====================================================

    def write_reference_haplotypes_merged(self, out_fa: str):
        """
        把 75% reference 样本的 hap 全部写到一个大 FASTA 里：
        data/ref_haps.fa

        >sampleA_hap1
        ...
        >sampleA_hap2
        ...
        """
        ref_samples = self.reference_samples
        out_path = Path(out_fa)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w") as f_out:
            for s in ref_samples:
                hap1, hap2 = self._build_sample_haplotypes(s)
                seq1 = "".join(hap1)
                seq2 = "".join(hap2)

                f_out.write(f">{s}_hap1\n{seq1}\n")
                f_out.write(f">{s}_hap2\n{seq2}\n")

                print(f"[INFO] appended haplotypes of {s} to {out_fa}")

        print(f"[OK] reference haplotypes merged -> {out_fa}")

    # =====================================================
    #  4b. 写 25% test：每个样本一个目录 + fasta
    # =====================================================

    def write_test_haplotypes_per_sample(self, line_width: int = 100):
        """
        为 25% test 样本写单独的 FASTA：
        data/reads/<sample>/<sample>.fa
        """
        for s in self.test_samples:
            hap1, hap2 = self._build_sample_haplotypes(s)

            sample_dir = self.read_root / s    # 比如 data/reads/tsk_0
            fasta_path = sample_dir / f"{s}.fa"

            self._write_two_haps_to_fasta(
                sample=s,
                hap1=hap1,
                hap2=hap2,
                out_fa=fasta_path,
                line_width=line_width,
            )

    # =====================================================
    #  5. 对 25% test 样本跑 wgsim（双端）
    # =====================================================

    def run_wgsim_for_tests(
        self,
        depth: float,
        read_len: int = 120,
        insert_mean: int = 300,
        insert_sd: int = 50,
        err: float = 0.001,
        seed: int = 42,
        wgsim_path: str = "wgsim",
    ):
        """
        只对 25% test 样本调用 wgsim，生成双端 reads：

        data/reads/<sample>/<sample>_R1.fq
        data/reads/<sample>/<sample>_R2.fq

        depth: 每个样本目标测序深度（近似 2L 上的平均深度）
        """
        L = self.L

        for s in self.test_samples:
            sample_dir = self.read_root / s
            fasta_path = sample_dir / f"{s}.fa"
            r1_path = sample_dir / f"{s}_R1.fq"
            r2_path = sample_dir / f"{s}_R2.fq"

            if not fasta_path.exists():
                print(f"[WARN] {fasta_path} not found, skip {s}")
                continue

            # 二倍体总长度 ~ 2L，每对 PE read 覆盖 ~2*read_len
            # 目标覆盖 depth * 2L
            # N_pairs ≈ depth * 2L / (2*read_len) = depth * L / read_len
            N_pairs = math.ceil(depth * L / read_len)
            print(f"[INFO] wgsim for {s}: depth={depth}, L={L}, N_pairs={N_pairs}")

            cmd = [
                wgsim_path,
                "-N", str(N_pairs),
                "-1", str(read_len),
                "-2", str(read_len),
                "-e", str(err),
                "-r", "0",
                "-R", "0",
                "-X", "0",
                "-d", str(insert_mean),
                "-s", str(insert_sd),
                "-S", str(seed),
                str(fasta_path),
                str(r1_path),
                str(r2_path),
            ]
            subprocess.run(cmd, check=True)


import pandas as pd
from typing import List, Tuple

def stratified_split_samples_pd(
    meta_path: str,
    ref_out: str,
    test_out: str,
    ref_ratio: float = 0.75,
    seed: int = 42,
    sample_col: str = "sample",
    pop_col: str = "population",
    sep: str = "\t",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    用 pandas 按群体(pop_col)分层，把样本划分为参考集(ref)和测试集(test)，
    并把样本及其对应种群一起写到输出文件中。

    输出文件格式（TSV）：
    sample_id <tab> population

    返回
    ----
    (ref_df, test_df) 两个 DataFrame，只包含 sample_col 和 pop_col 两列
    """
    # 读入样本表
    df = pd.read_csv(meta_path, sep=sep)

    # 为了可复现，先整体打乱
    df = df.sample(frac=1.0, random_state=seed).reset_index(drop=True)

    ref_indices = []
    test_indices = []

    # 按群体分组，分别做 75/25 划分
    for pop, sub in df.groupby(pop_col):
        n = len(sub)
        n_ref = int(round(n * ref_ratio))

        ref_indices.extend(sub.index[:n_ref].tolist())
        test_indices.extend(sub.index[n_ref:].tolist())

        print(f"[INFO] pop={pop}: total={n}, ref={n_ref}, test={n - n_ref}")

    ref_df = df.loc[ref_indices, [sample_col, pop_col]].reset_index(drop=True)
    test_df = df.loc[test_indices, [sample_col, pop_col]].reset_index(drop=True)

    # 写出：两列 sample_id + population
    ref_df.to_csv(ref_out, index=False)
    test_df.to_csv(test_out, index=False)

    print(f"[OK] reference list -> {ref_out}")
    print(f"[OK] test list      -> {test_out}")

    return ref_df, test_df
