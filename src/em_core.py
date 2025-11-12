# -*- coding: utf-8 -*-
import math
import numpy as np
import pandas as pd

class EMEstimator(object):
    """
    只负责：
      - 读取单样本 read-hap 配对 TSV
      - 预计算 log q_ij
      - 运行 EM 得到 p_i (haplotype mixing proportions)
    """

    def __init__(self, K, alpha=0.005, max_iter=300, tol=1e-5, verbose=True):
        self.K = int(K)
        self.alpha = float(alpha)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.verbose = bool(verbose)

    def load_alignment_tsv(self, tsv_path):
        """
        期望列：read_id, hap_index, d_ij, n_j
        允许多行同一 read 对应多个 hap（–all 比对）
        """
        df = pd.read_csv(tsv_path, sep=r"[,\t]", engine="python")
        need = ["read_id", "hap_index", "d_ij", "n_j"]
        for c in need:
            if c not in df.columns:
                raise ValueError("TSV missing column: %s" % c)

        df = df.copy()
        df["hap_index"] = df["hap_index"].astype(int)
        df["d_ij"] = df["d_ij"].astype(int)
        df["n_j"] = df["n_j"].astype(int)

        # 为每条 read 指派递增索引
        df["read_idx"] = pd.factorize(df["read_id"])[0].astype(int)
        return df

    def precompute_log_q(self, df):
        """
        q_ij = alpha^{d_ij} * (1-alpha)^{n_j - d_ij}
        返回与 df 对齐的 log_q 向量
        """
        d = df["d_ij"].to_numpy(dtype=np.int64, copy=False)
        n = df["n_j"].to_numpy(dtype=np.int64, copy=False)

        a = min(max(self.alpha, 1e-12), 1.0 - 1e-12)
        loga = math.log(a)
        log1a = math.log(1.0 - a)

        log_q = d * loga + (n - d) * log1a
        return log_q.astype(np.float64, copy=False)

    def run_em(self, df, log_q):
        """
        输入：
          df: load_alignment_tsv 的返回
          log_q: precompute_log_q 的返回
        输出：
          p:  np.ndarray[K]  （非负、和为 1）
          ll: np.ndarray[<=max_iter]（每次迭代的对数似然）
        """
        K = self.K

        # 将记录按 read_idx 排序，便于按 read 切片遍历
        order = np.argsort(df["read_idx"].to_numpy())
        df_sorted = df.iloc[order].reset_index(drop=True)
        log_q_sorted = log_q[order]

        read_ids = df_sorted["read_idx"].to_numpy()
        bounds = np.flatnonzero(np.r_[True, read_ids[1:] != read_ids[:-1], True])

        # p 初始为均匀
        p = np.full(K, 1.0 / K, dtype=np.float64)
        ll_hist = []

        for it in range(1, self.max_iter + 1):
            numer = np.zeros(K, dtype=np.float64)  # \sum_j r_ij
            ll_val = 0.0

            for bi in range(len(bounds) - 1):
                s = bounds[bi]
                e = bounds[bi + 1]

                hap_idx = df_sorted["hap_index"].to_numpy()[s:e]
                lq = log_q_sorted[s:e]

                # log p_i + log q_ij
                log_pi = np.log(p[hap_idx])
                log_terms = log_pi + lq

                # log-sum-exp
                m = np.max(log_terms)
                denom = m + math.log(np.sum(np.exp(log_terms - m)))
                ll_val += denom

                # r_ij = exp(log_terms - denom)
                rij = np.exp(log_terms - denom)
                np.add.at(numer, hap_idx, rij)

            total = numer.sum()
            if total <= 0.0 or not np.isfinite(total):
                new_p = np.full(K, 1.0 / K, dtype=np.float64)
            else:
                new_p = numer / total

            delta = np.max(np.abs(new_p - p))
            p = new_p
            ll_hist.append(ll_val)

            if self.verbose and (it % 10 == 0 or it == 1):
                print("[EM] iter=%d  ll=%.6f  delta=%.3e" % (it, ll_val, delta))

            if delta < self.tol:
                if self.verbose:
                    print("[EM] converged at iter=%d  ll=%.6f" % (it, ll_val))
                break

        # 归一化保护
        p = p.clip(min=0.0)
        s = p.sum()
        if s <= 0:
            p[:] = 1.0 / K
        else:
            p /= s

        return p, np.asarray(ll_hist, dtype=np.float64)
