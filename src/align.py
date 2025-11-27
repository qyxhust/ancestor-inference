import subprocess
from pathlib import Path

import pandas as pd
import pysam


class EMAlignment(object):
    """
    功能：
    1) 用 75% 参考 hap FASTA 建 bowtie2 索引
    2) 对 25% test 样本的 reads (R1/R2) 对齐到参考 hap
    3) 从 BAM 中抽取 (sample, read, hap, mismatch, read_len) 写成 EM 可用的 TSV
    """

    def __init__(
        self,
        ref_fasta_path,
        test_samples_csv,
        reads_root,
        out_root,
        index_prefix=None,
    ):
        """
        ref_fasta_path : 参考 hap 的 multi-FASTA，例如 data/ref_haps/ref_haps.fa
                        每条记录是一条 haplotype
        test_samples_csv : 包含一列 'sample' 的 CSV（你之前分好的 25% 样本列表）
        reads_root      : 单样本 reads 根目录：
                            reads_root/<sample>/<sample>_R1.fq
                            reads_root/<sample>/<sample>_R2.fq
        out_root        : 输出根目录：
                            out_root/align/     存 BAM
        index_prefix    : bowtie2 索引前缀（默认= ref_fasta 去掉后缀再加 _index）
        """
        self.ref_fasta_path = Path(ref_fasta_path)
        self.test_samples_csv = Path(test_samples_csv)
        self.reads_root = Path(reads_root)

        self.out_root = Path(out_root)
        self.align_dir = self.out_root / "bam"
        self.out_root.mkdir(parents=True, exist_ok=True)
        self.align_dir.mkdir(parents=True, exist_ok=True)

        if index_prefix is None:
            # e.g. data/ref_haps/ref_haps.fa -> data/ref_haps/ref_haps_index
            self.index_prefix = str(self.ref_fasta_path).rsplit(".", 1)[0] + "_index"
        else:
            self.index_prefix = str(index_prefix)

        # 读 test 样本列表
        df = pd.read_csv(self.test_samples_csv)
        if "sample" not in df.columns:
            raise ValueError("test_samples_csv 必须包含列 'sample'")
        self.test_samples = df["sample"].tolist()

    # ---------- 1. 建 bowtie2 索引 + hap 索引表 ----------

    def build_bowtie2_index(self, bowtie2_build="bowtie2-build"):
        """
        用参考 hap FASTA 建 bowtie2 索引。
        如果已存在 *.1.bt2，就跳过。
        """
        bt2_file = self.index_prefix + ".1.bt2"
        if Path(bt2_file).exists():
            print("[INFO] bowtie2 index 已存在，跳过构建：", self.index_prefix)
            return

        cmd = [bowtie2_build, str(self.ref_fasta_path), self.index_prefix]
        print("[INFO] Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        print("[OK] bowtie2 index built:", self.index_prefix)

    def write_hap_index_table(self, out_csv=None):
        """
        从 ref_fasta 的顺序写出：
        hap_index,hap_id
        """
        if out_csv is None:
            out_csv = self.out_root / "ref_haps_index.csv"
        else:
            out_csv = Path(out_csv)

        fasta = pysam.FastaFile(str(self.ref_fasta_path))
        refs = list(fasta.references)
        rows = [{"hap_index": i, "hap_id": name} for i, name in enumerate(refs)]
        fasta.close()

        df = pd.DataFrame(rows)
        df.to_csv(out_csv, index=False)
        print("[OK] hap index table ->", out_csv)
        return df

    # ---------- 2. 对样本做 bowtie2 对齐 ----------

    def align_one_sample(
        self,
        sample_id,
        bowtie2_path="bowtie2",
        samtools_path="samtools",
        threads=4,
    ):
        """
        对单个 test 样本的 R1/R2 做 bowtie2 --all 对齐，
        输出 sorted BAM： out_root/align/<sample>.sorted.bam
        """
        r1 = self.reads_root / sample_id / f"{sample_id}_R1.fq"
        r2 = self.reads_root / sample_id / f"{sample_id}_R2.fq"

        if not r1.exists() or not r2.exists():
            raise FileNotFoundError(f"找不到 fastq：{r1} 或 {r2}")

        sam_path = self.align_dir / f"{sample_id}.sam"
        bam_path = self.align_dir / f"{sample_id}.sorted.bam"

        # 1) bowtie2 -> SAM
        cmd_bt = [
            bowtie2_path,
            "--all",
            "--very-sensitive",
            "-p",
            str(threads),
            "-x",
            self.index_prefix,
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-S",
            str(sam_path),
        ]
        print("[INFO] bowtie2 for sample:", sample_id)
        print("       " + " ".join(cmd_bt))
        subprocess.run(cmd_bt, check=True)

        # 2) SAM -> sorted BAM
        cmd_view = [samtools_path, "view", "-bS", str(sam_path)]
        cmd_sort = [samtools_path, "sort", "-o", str(bam_path)]
        print("[INFO] samtools view|sort for sample:", sample_id)
        p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd_sort, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
        if p2.returncode != 0:
            raise RuntimeError(f"samtools sort 失败 for sample {sample_id}")

        # 可选：删掉 SAM
        try:
            sam_path.unlink()
        except OSError:
            pass

        print("[OK] BAM written:", bam_path)
        return bam_path

    def align_all_samples(self, bowtie2_path="bowtie2", samtools_path="samtools", threads=4):
        """
        对所有 test 样本批量对齐。
        返回 {sample_id: bam_path}。
        """
        bam_paths = {}
        for s in self.test_samples:
            bam_paths[s] = self.align_one_sample(
                s,
                bowtie2_path=bowtie2_path,
                samtools_path=samtools_path,
                threads=threads,
            )
        return bam_paths