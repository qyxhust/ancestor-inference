
from pathlib import Path
import pysam
import numpy as np
import subprocess
import pysam
import numpy as np
from pathlib import Path
import textwrap


class SequencingRead:
    def __init__(self, vcf_path: str, L: int, fasta_path: str, chrom: str = None):
        """
        vcf_path: msprime 生成的 VCF
        L: 序列长度（需要 >= VCF 中最大坐标）
        out_dir: 输出目录（FASTA / 以后 wgsim 输出）
        chrom: 要处理的染色体名；如果 None，就用 VCF 中第一条
        """
        self.vcf_path = vcf_path
        self.L = L
        self.fasta_path = Path(fasta_path)


        self.vcf = pysam.VariantFile(vcf_path)
        self.samples = list(self.vcf.header.samples)

        self.chrom = "1"

    def _build_ref_array(self):
        """
        构建参考序列数组：
        1. 初始全 A
        2. 在所有变异位点上用 REF 覆盖
        返回: numpy array, 长度 L, 元素是单个字符 'A','C','G','T'
        """
        ref_array = np.full(self.L, "A", dtype="U1")

        # 重新打开 VCF，避免 self.vcf 已经被迭代过的问题
        vcf = pysam.VariantFile(self.vcf_path)
        for rec in vcf.fetch(self.chrom):
            # 这里只处理 SNP（单碱基 REF）
            if len(rec.ref) != 1:
                continue
            pos = rec.pos - 1  # VCF 1-based -> 0-based
            if 0 <= pos < self.L:
                ref_array[pos] = rec.ref

        return ref_array
    
    def write_all_sample_haplotypes_fasta(self,
                                        line_width: int = 100):
        """
        基于 VCF 和参考序列，为每个样本构建两条 haplotype，
        不在类里保存大数组，按“一个样本一处理一写入”的方式写入同一个 FASTA：

        >sample_hap1
        SEQUENCE...
        >sample_hap2
        SEQUENCE...
        """
        # 1. 构建参考序列数组（全 A -> 用 REF 覆盖 SNP 位点）
        base_ref_array = self._build_ref_array()
        vcf = pysam.VariantFile(self.vcf_path)

        with open(self.fasta_path, "w") as f_out:

            # ===== 外层：按样本遍历 =====
            for s in self.samples:
                # 为当前样本准备两条 haplotype（初始都是 REF）
                hap1 = base_ref_array.copy()
                hap2 = base_ref_array.copy()

                # ===== 内层：遍历所有 SNP 位点，更新当前样本的 hap1/2 =====
                for rec in vcf.fetch(self.chrom):
                    if rec.contig != self.chrom:
                        continue
                    # 只处理 SNP
                    if len(rec.ref) != 1:
                        continue

                    pos = rec.pos - 1
                    if not (0 <= pos < self.L):
                        continue

                    alts = rec.alts  # 例如 ('C',) 或 ('C','G')
                    if not alts:
                        continue

                    sample_data = rec.samples[s]
                    gt = sample_data.get("GT", None)  # 例如 (0,1) / (1,1) / (0,0)
                    if gt is None:
                        continue

                    # hap_idx: 0 -> hap1, 1 -> hap2
                    for hap_idx, allele_index in enumerate(gt):
                        if allele_index is None or allele_index == 0:
                            # None 或 0（REF）都不用改
                            continue

                        # ALT 等位基因（1 -> alts[0], 2 -> alts[1], ...）
                        if allele_index - 1 >= len(alts):
                            continue
                        alt_base = alts[allele_index - 1]
                        if len(alt_base) != 1:
                            # 非 SNP（indel）这里直接跳过
                            continue

                        if hap_idx == 0:
                            hap1[pos] = alt_base
                        else:  # hap_idx == 1
                            hap2[pos] = alt_base

                # ===== 当前样本的所有 SNP 都更新完，写入 FASTA =====
                seq1 = "".join(hap1)
                seq2 = "".join(hap2)

                f_out.write(f">{s}_hap1\n")
                for chunk in textwrap.wrap(seq1, line_width):
                    f_out.write(chunk + "\n")

                f_out.write(f">{s}_hap2\n")
                for chunk in textwrap.wrap(seq2, line_width):
                    f_out.write(chunk + "\n")

        print(f"[INFO] Written all samples haplotypes to {self.fasta_path}")

    
    #     使用 subprocess 调用 wgsim，对 all_samples.fa 生成 fastq
    def run_wgsim(self,
                fasta_name: str = "all_samples.fa",
                out_prefix: str = "reads",
                N: int = 100000,
                read_len: int = 150,
                err: float = 0.001,
                seed: int = 42,
                wgsim_path: str = "wgsim"):
        """
        调用 wgsim 对 FASTA 做测序模拟（成对 reads）。

        参数：
        - fasta_name: 之前写出的 FASTA 文件名（在 out_dir 下面）
        - out_prefix: 输出 fastq 的前缀，生成 out_prefix_R1.fq / out_prefix_R2.fq
        - N: 需要模拟的 read 对数（wgsim -N）
        - read_len: read 长度（-1 和 -2 都用这个长度）
        - err: 测序错误率 (-e)
        - seed: 随机数种子 (-S)
        - wgsim_path: wgsim 可执行文件路径（在 PATH 里就直接用 "wgsim"）
        """
        fasta_path = self.out_dir / fasta_name
        r1_path = self.out_dir / f"{out_prefix}_R1.fq"
        r2_path = self.out_dir / f"{out_prefix}_R2.fq"

        cmd = [
            wgsim_path,
            "-N", str(N),
            "-1", str(read_len),
            "-2", str(read_len),
            "-e", str(err),
            "-S", str(seed),
            str(fasta_path),
            str(r1_path),
            str(r2_path),
        ]

        print("[INFO] Running wgsim:")
        print("      " + " ".join(cmd))

        subprocess.run(cmd, check=True)

        print(f"[OK] wgsim finished. Output:")
        print(f"     {r1_path}")
        print(f"     {r2_path}")
