import numpy as np
import pandas as pd
from pathlib import Path


def generate_pool_weights_csv(
    classes,
    out_csv: str,
    n_pools=None,
    pool_ids=None,
    seed: int = 42,
):
    """
    生成多个 pool 的类权重（每行一个 pool），写到 out_csv。

    用法：
    - 只给 n_pools：自动生成 Pool1..PoolN
    - 或者直接给 pool_ids 列表
    """
    if pool_ids is None:
        if n_pools is None:
            raise ValueError("需要提供 n_pools 或 pool_ids 其中之一")
        pool_ids = [f"Pool{i}" for i in range(1, n_pools + 1)]

    rng = np.random.default_rng(seed)
    rows = []

    for pool_id in pool_ids:
        raw = rng.random(len(classes))
        weights = raw / raw.sum()  # 每个 pool 内权重和为 1

        row = {"pool_id": pool_id}
        for cls, w in zip(classes, weights):
            row[cls] = w
        rows.append(row)

    df = pd.DataFrame(rows, columns=["pool_id"] + list(classes))
    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[OK] pool weights -> {out_path}")
    return 

import numpy as np
import pandas as pd
from pathlib import Path


class Pooling:
    """
    只负责：给定单个 pool 的权重 + 样本信息 + 单样本 reads，
    生成这个 pool 的混合双端测序文件。
    """

    def __init__(self, pools_dir="data/pools"):
        self.pools_dir = Path(pools_dir)
        self.pools_dir.mkdir(parents=True, exist_ok=True)

    # ---- 工具：从一个样本的 R1/R2 中抽 n_pairs 对 reads ----
    @staticmethod
    def _sample_pairs_from_fastq(
        r1_path,
        r2_path,
        n_pairs,
        out_r1,
        out_r2,
        rng,
    ):
        """
        从一个样本的 R1/R2 中随机抽 n_pairs 对 read，写到 out_r1/out_r2。
        简单实现：两趟扫描。
        """
        r1_path = Path(r1_path)
        r2_path = Path(r2_path)

        # 第一次：统计总 pair 数
        with r1_path.open() as f1:
            total_pairs = sum(1 for _ in f1) // 4

        if total_pairs == 0 or n_pairs <= 0:
            return

        if n_pairs >= total_pairs:
            # 不下采样，直接全部拷贝
            with r1_path.open() as f1, r2_path.open() as f2:
                for line in f1:
                    out_r1.write(line)
                for line in f2:
                    out_r2.write(line)
            return

        # 随机选择要保留的 pair index
        chosen_idx = rng.choice(total_pairs, size=n_pairs, replace=False)
        chosen_idx.sort()
        chosen_set = set(chosen_idx.tolist())

        # 第二次：真正写出选中的 pair
        with r1_path.open() as f1, r2_path.open() as f2:
            pair_idx = 0
            while True:
                r1_block = [f1.readline() for _ in range(4)]
                r2_block = [f2.readline() for _ in range(4)]

                if r1_block[0] == "" or r2_block[0] == "":
                    break

                if pair_idx in chosen_set:
                    for line in r1_block:
                        out_r1.write(line)
                    for line in r2_block:
                        out_r2.write(line)

                pair_idx += 1

    # ---- 核心：根据“这一行权重 + 样本信息 + depth”生成单个 pool ----
    def create_single_pool(
        self,
        pool_id,
        classes,
        class_weights,        # dict 或 Series：键为 classes，值为权重
        sample_meta_csv,
        reads_root,
        L,
        depth,
        hap_per_pool,
        sigma_frac=0.2,
        read_len=150,
        seed=123,
    ):
        """
        pool_id       : 当前 pool 名称（用于输出文件名）
        classes       : ['YRI','CEU','NEA']（列的顺序）
        class_weights : 这一行的权重，比如 {'YRI':0.5,'CEU':0.3,'NEA':0.2}
        sample_meta_csv: test 样本信息 CSV（包含 'sample','population'）
        reads_root    : 单样本 reads 根目录，reads_root/<sample>/<sample>_R1.fq/_R2.fq
        L             : 单倍体长度
        depth         : 全局测序深度（与 wgsim 一致）
        hap_per_pool  : 每个 pool 单倍体总数（偶数），样本数 = hap_per_pool/2
        sigma_frac    : per-sample read 数的标准差比例，σ = sigma_frac * μ
        read_len      : read 长度
        seed          : 随机种子
        """
        rng = np.random.default_rng(seed)

        # 1. 读样本信息
        meta = pd.read_csv(sample_meta_csv)
        if "sample" not in meta.columns or "population" not in meta.columns:
            raise ValueError("sample_meta_csv 需要包含列 'sample','population'")

        reads_root = Path(reads_root)

        # 2. pool 中样本数
        if hap_per_pool % 2 != 0:
            raise ValueError("hap_per_pool 必须是偶数（每个样本两条 hap）")
        pool_sample_count = hap_per_pool // 2

        # 3. 从 class_weights 里按 classes 提取权重，并决定每类多少个样本
        weights = np.array([float(class_weights[cls]) for cls in classes])
        weights = weights / weights.sum()
        n_raw = weights * pool_sample_count
        n_cls = np.floor(n_raw).astype(int)
        diff = pool_sample_count - n_cls.sum()
        if diff > 0:
            order = np.argsort(-weights)
            for idx in order[:diff]:
                n_cls[idx] += 1
        class_to_n = {cls: int(n) for cls, n in zip(classes, n_cls)}
        print("[INFO]", pool_id, "pool_sample_count=", pool_sample_count,
              "class_to_n=", class_to_n)

        # 4. 在每个类中随机挑样本
        selected_samples = []
        for cls in classes:
            need = class_to_n[cls]
            if need <= 0:
                continue
            sub = meta[meta["population"] == cls]
            if len(sub) < need:
                raise ValueError(
                    "%s: class %s 需要 %d 个样本，但只有 %d 个"
                    % (pool_id, cls, need, len(sub))
                )
            chosen = sub.sample(
                n=need,
                replace=False,
                random_state=rng.integers(0, 2**32 - 1),
            )
            selected_samples.extend(chosen["sample"].tolist())

        n_selected = len(selected_samples)
        if n_selected == 0:
            raise ValueError("%s: 没选出任何样本" % pool_id)
        print("[INFO]", pool_id, "selected", n_selected, "samples:", selected_samples)

        # 5. pool 总 read 对数：沿用全局 depth
        N_pool = int(np.ceil(depth * L / read_len))
        mu = float(N_pool) / n_selected
        sigma = sigma_frac * mu
        print("[INFO]", pool_id, "N_pool≈%d, μ≈%.1f, σ≈%.1f" % (N_pool, mu, sigma))

        # 6. 为每个样本抽 Ns（正态），截断到 [0, 该样本的总 pair 数]
        Ns = {}
        total_actual = 0
        for s in selected_samples:
            r1_path = reads_root / s / ("%s_R1.fq" % s)
            if not r1_path.exists():
                raise FileNotFoundError("%s: %s 不存在" % (pool_id, str(r1_path)))
            with r1_path.open() as f1:
                n_total = sum(1 for _ in f1) // 4

            val = rng.normal(loc=mu, scale=sigma)
            N_s = int(round(val))
            if N_s < 0:
                N_s = 0
            if N_s > n_total:
                N_s = n_total

            Ns[s] = N_s
            total_actual += N_s

        print("[INFO]", pool_id,
              "total Ns after normal sampling =", total_actual)

        # 7. 抽样并写出 pool R1/R2
        out_r1_path = self.pools_dir / ("%s_R1.fq" % pool_id)
        out_r2_path = self.pools_dir / ("%s_R2.fq" % pool_id)

        with out_r1_path.open("w") as out_r1, out_r2_path.open("w") as out_r2:
            for s in selected_samples:
                N_s = Ns[s]
                if N_s <= 0:
                    continue

                r1_path = reads_root / s / ("%s_R1.fq" % s)
                r2_path = reads_root / s / ("%s_R2.fq" % s)
                if (not r1_path.exists()) or (not r2_path.exists()):
                    print("[WARN]", pool_id, "sample", s,
                          "fastq 缺失，跳过")
                    continue

                print("[INFO]", pool_id, "sample =", s, ", N_pairs =", N_s)
                sub_rng = np.random.default_rng(rng.integers(0, 2**32 - 1))
                self._sample_pairs_from_fastq(
                    r1_path=r1_path,
                    r2_path=r2_path,
                    n_pairs=N_s,
                    out_r1=out_r1,
                    out_r2=out_r2,
                    rng=sub_rng,
                )

        print("[OK]", pool_id, "pool fastq written:",
              str(out_r1_path), str(out_r2_path))


        