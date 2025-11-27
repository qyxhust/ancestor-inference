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
from src.pool import generate_pool_weights_csv  # 负责：生成 pool 权重 CSV
from scripts.run_pool import run_pool          # 你的函数，内部已内置 cfg_path
from scripts.run_align import run_align        # 你的函数，内部已内置 cfg_path
from scripts.run_em import run_em_for_all_tests             # 你的函数，内部已内置 cfg_path



def main():
    # # 1) 运行模拟，生成 VCF + 标签
    # run_simulate(model_type="hard")
    # print("[OK] Pipeline complete.")

    # # 2) 分层划分样本
    # stratified_split_samples_pd(
    #     meta_path="/space/s1/qyx/data/simulate/labels.tsv",
    #     ref_out="/space/s1/qyx/data/ref_samples.csv",
    #     test_out="/space/s1/qyx/data/test_samples.csv",
    #     ref_ratio=0.75,
    #     seed=42,
    #     sep="\t",
    # )

    # run_read()
    # print("[OK] Read simulation complete.")

    # # 3) 生成 pool 权重 CSV
    # classes = ["YRI", "CEU", "Neandertal"]
    # pool_weight_csv = "data/pools/pool_weights.csv"

    # weights_df = generate_pool_weights_csv(
    #     classes=classes,
    #     out_csv=pool_weight_csv,
    #     n_pools=10,      # 👈 没给 pool_ids，就自动用 Pool1..Pool4
    #     seed=42,
    # )

    # # 4) 运行 pool 混样，生成混合测序文件
    # run_pool()
    # print("[OK] Pooling complete.")

    # 5) 运行对齐，生成 EM 输入文件
    run_align()
    # print("[OK] Alignment complete.")

    # 6) 运行 EM，生成最终结果
    # 你当前工程默认结构（可按需修改）
    # tsv_dir        = "data/em"
    # hap_index_csv  = "data/align/ref_haps_index.csv"
    # ref_meta_csv   = "data/ref_samples.csv"
    # test_list_csv  = "data/test_samples.csv"
    # out_csv        = "data/em/em_results.csv"

    # # 跑全部测试样本，逐样本写结果（追加）
    # run_em_for_all_tests(
    #     tsv_dir=tsv_dir,
    #     test_list_csv=test_list_csv,
    #     hap_index_csv=hap_index_csv,
    #     ref_meta_csv=ref_meta_csv,
    #     out_csv=out_csv,
    #     alpha=0.005,
    #     max_iter=100,
    #     tol=1e-5,
    #     verbose=True
    # )
    # print("[PIPELINE] EM 全部完成 -> %s" % out_csv)

if __name__ == "__main__":
    main()
