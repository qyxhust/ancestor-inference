from src.align import EMAlignment

def run_align():
    # 75% 参考 hap 的 multi-FASTA
    ref_fasta = "data/ref_haps.fa"

    # 25% test 样本列表（有一列 sample）
    test_samples_csv = "data/test_samples.csv"

    # 单样本 reads 根目录：<reads_root>/<sample>/<sample>_R1.fq/_R2.fq
    reads_root = "data/reads"

    # 对齐 & EM 输入的输出目录
    out_root = "data/align"

    aligner = EMAlignment(
        ref_fasta_path=ref_fasta,
        test_samples_csv=test_samples_csv,
        reads_root=reads_root,
        out_root=out_root,
    )

    # 1. 建 bowtie2 索引 + hap_index 表
    aligner.build_bowtie2_index(bowtie2_build="bowtie2-build")
    hap_index_df = aligner.write_hap_index_table()  # data/em/ref_haps_index.csv

    # 2. 对所有 test 样本做对齐，生成 sorted BAM
    aligner.align_all_samples(
        bowtie2_path="bowtie2",
        samtools_path="samtools",
        threads=8,
    )

    # 3. 从 BAM 生成 EM 需要的 read–hap mismatch TSV
    aligner.em_input_all_samples(hap_index_df=hap_index_df)
