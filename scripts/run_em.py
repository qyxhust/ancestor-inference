# -*- coding: utf-8 -*-
import os
from pathlib import Path
import pandas as pd
import numpy as np

from src.em_core import EMEstimator

# ---------- 基础工具 ----------
def _ensure_parent_dir(path_str):
    p = Path(path_str)
    p.parent.mkdir(parents=True, exist_ok=True)

def _populations_from_meta(ref_meta_csv):
    meta = pd.read_csv(ref_meta_csv)
    if "population" not in meta.columns:
        raise ValueError("ref_meta_csv must contain column 'population'")
    pops = list(pd.unique(meta["population"]))
    pops.sort()
    return pops

def _aggregate_to_populations(hap_prob_df, ref_meta_csv):
    """
    输入 hap_prob_df: ['sample','hap_index','p','hap_id']
    其中 hap_id 形如 "tsk_123_hap1"
    ref_meta: 两列 ['sample','population']
    返回 long 格式：['sample','population','p_population']
    """
    meta = pd.read_csv(ref_meta_csv)
    if "sample" not in meta.columns or "population" not in meta.columns:
        raise ValueError("ref_meta_csv must contain columns: sample, population")

    tmp = hap_prob_df.copy()
    tmp["sample_ref"] = tmp["hap_id"].str.replace(r"_hap[12]$", "", regex=True)

    tmp = tmp.merge(meta, left_on="sample_ref", right_on="sample", how="left", suffixes=("", "_ref"))
    if tmp["population"].isna().any():
        missing = tmp.loc[tmp["population"].isna(), "sample_ref"].unique()
        raise ValueError("Some reference haplotypes have no population label: %s" % list(missing))

    grp = tmp.groupby(["sample", "population"], as_index=False)["p"].sum()
    grp = grp.rename(columns={"p": "p_population"})
    return grp

def _aggregate_to_populations_wide(hap_prob_df, ref_meta_csv, pop_order=None):
    long_df = _aggregate_to_populations(hap_prob_df, ref_meta_csv)
    wide = long_df.pivot(index="sample", columns="population", values="p_population").fillna(0.0)
    wide = wide.reset_index().rename_axis(None, axis=1)

    if pop_order is None:
        pop_order = _populations_from_meta(ref_meta_csv)

    # 保证所有群体列都在，且顺序固定
    for c in pop_order:
        if c not in wide.columns:
            wide[c] = 0.0

    cols = ["sample"] + pop_order
    return wide[cols]

# ---------- 单样本：跑 EM -> 聚合到种群 -> 返回一行 ----------
def run_em_file_to_pop_row(tsv_path, hap_index_csv, ref_meta_csv,
                           alpha=0.005, max_iter=300, tol=1e-5, verbose=True, pop_order=None):
    """
    返回 DataFrame 一行： sample,<POP1>,<POP2>,<POP3>,...
    """
    hap_df = pd.read_csv(hap_index_csv)
    if "hap_index" not in hap_df.columns or "hap_id" not in hap_df.columns:
        raise ValueError("hap_index_csv must contain columns: hap_index,hap_id")
    K = int(hap_df["hap_index"].max()) + 1

    est = EMEstimator(K=K, alpha=alpha, max_iter=max_iter, tol=tol, verbose=verbose)
    df = est.load_alignment_tsv(tsv_path)
    log_q = est.precompute_log_q(df)
    p, _ = est.run_em(df, log_q)

    sample_name = Path(tsv_path).name.replace("_read_hap.tsv", "")
    hap_prob_df = pd.DataFrame({
        "sample": sample_name,
        "hap_index": np.arange(K, dtype=int),
        "p": p
    }).merge(hap_df, on="hap_index", how="left")

    wide = _aggregate_to_populations_wide(hap_prob_df, ref_meta_csv, pop_order=pop_order)
    return wide  # 一行 DataFrame

# ---------- 批量跑：每个样本跑完立刻落盘 ----------
def run_em_for_all_tests(tsv_dir, test_list_csv, hap_index_csv, ref_meta_csv, out_csv,
                         alpha=0.005, max_iter=300, tol=1e-5, verbose=True):
    """
    - tsv_dir:            data/align/tsv/
    - test_list_csv:      data/test_samples.csv（两列：sample,population）
    - hap_index_csv:      data/align/ref_haps_index.csv（两列：hap_index,hap_id）
    - ref_meta_csv:       data/ref_samples.csv（两列：sample,population）
    - out_csv:            data/em/em_results.csv（每行：sample,<POP1>,<POP2>,...）
    """
    tsv_dir = Path(tsv_dir)
    out_csv = Path(out_csv)
    _ensure_parent_dir(out_csv)

    tests = pd.read_csv(test_list_csv)
    if "sample" not in tests.columns:
        raise ValueError("test_list_csv must contain column 'sample'")

    pop_order = _populations_from_meta(ref_meta_csv)

    # 如果文件已存在，后续只追加不存在的样本；否则写表头
    header_needed = not out_csv.exists()
    if header_needed:
        # 先写一个空表头
        hdr_df = pd.DataFrame(columns=["sample"] + pop_order)
        hdr_df.to_csv(out_csv, index=False)

    # 已有结果的样本集合
    done = set()
    if out_csv.exists():
        try:
            exist = pd.read_csv(out_csv)
            if "sample" in exist.columns:
                done = set(exist["sample"].astype(str).tolist())
        except Exception:
            done = set()

    total = len(tests)
    for i, s in enumerate(tests["sample"].astype(str).tolist(), 1):
        tsv_path = tsv_dir / ("%s_read_hap.tsv" % s)
        if not tsv_path.exists():
            print("[EM-BATCH] (%d/%d) %s  对齐文件不存在，跳过" % (i, total, s))
            continue

        if s in done:
            print("[EM-BATCH] (%d/%d) %s  已存在结果，跳过" % (i, total, s))
            continue

        print("[EM-BATCH] (%d/%d) %s  正在运行..." % (i, total, s))
        try:
            wide_row = run_em_file_to_pop_row(
                str(tsv_path), str(hap_index_csv), str(ref_meta_csv),
                alpha=alpha, max_iter=max_iter, tol=tol, verbose=verbose, pop_order=pop_order
            )
            # 追加写入
            with open(out_csv, "a") as f:
                wide_row.to_csv(f, header=False, index=False)
            done.add(s)
            print("[EM-BATCH] (%d/%d) %s  完成并写入" % (i, total, s))
        except Exception as e:
            print("[EM-BATCH] (%d/%d) %s  失败: %s" % (i, total, s, str(e)))
