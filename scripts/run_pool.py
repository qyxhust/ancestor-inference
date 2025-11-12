import pandas as pd
from src.pool import Pooling
import yaml
from pathlib import Path

def run_pool():
    cfg = yaml.safe_load(Path("config/default.yaml").read_text())

    classes = ["YRI", "CEU", "Neandertal"]
    pool_weight_csv = "data/pools/pool_weights.csv"

    weights_df = pd.read_csv(pool_weight_csv)

    P = Pooling(pools_dir="data/pools")

    sample_meta_csv = "data/test_samples.csv"   # 有列 sample,population
    reads_root = "data/reads"             # 单样本 reads 根目录    
    hap_per_pool = 50                    # 每个 pool 30 hap -> 15 个样本

    for i, row in weights_df.iterrows():
        pool_id = row["pool_id"]                   # Pool1 / Pool2 / ...
        class_weights = {cls: row[cls] for cls in classes}

        P.create_single_pool(
            pool_id=pool_id,
            classes=classes,
            class_weights=class_weights,
            sample_meta_csv=sample_meta_csv,
            reads_root=reads_root,
            L=cfg["msprime"]["l"],
            depth=cfg["reads"]["depth"],
            hap_per_pool=hap_per_pool,
            sigma_frac=0.1,           # 波动程度可以自己调
            read_len=120,
            seed=42 + i,             # 每个 pool 一个不同 seed
        )







