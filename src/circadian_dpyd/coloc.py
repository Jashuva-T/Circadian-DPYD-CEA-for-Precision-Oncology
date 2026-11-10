"""Approximate-Bayes-factor colocalization (Giambartolomei et al. 2014, coloc.abf).

Needs two variant-level summary-statistic tables over the SAME variants (columns:
snp, beta, se). Raw GTEx / rhythmicity data are NOT bundled; see data/README.md.
The PP4 reported in the abstract (0.83) is a config input, not recomputed here.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.special import logsumexp

def log_abf(beta, se, w=0.15):
    """Wakefield log approximate Bayes factor, prior sd w on effect size."""
    z2 = (np.asarray(beta) / np.asarray(se)) ** 2
    v = np.asarray(se) ** 2
    r = w**2 / (w**2 + v)
    return 0.5 * (np.log(1 - r) + r * z2)

def coloc_abf(d1: pd.DataFrame, d2: pd.DataFrame, p1=1e-4, p2=1e-4, p12=1e-5, w1=0.15, w2=0.15) -> dict:
    m = d1.merge(d2, on="snp", suffixes=("_1", "_2"))
    if m.empty:
        raise ValueError("no overlapping variants")
    l1 = log_abf(m.beta_1, m.se_1, w1)
    l2 = log_abf(m.beta_2, m.se_2, w2)
    s1, s2 = logsumexp(l1), logsumexp(l2)
    l_h0 = 0.0
    l_h1 = np.log(p1) + s1
    l_h2 = np.log(p2) + s2
    # H3: distinct causal variants (i != j); H4: shared causal variant
    mat = l1[:, None] + l2[None, :]
    diag = np.diag(mat)
    off = mat[~np.eye(len(m), dtype=bool)]
    l_h3 = np.log(p1) + np.log(p2) + (logsumexp(off) if off.size else -np.inf)
    l_h4 = np.log(p12) + logsumexp(diag)
    ls = np.array([l_h0, l_h1, l_h2, l_h3, l_h4])
    pp = np.exp(ls - logsumexp(ls))
    return {f"PP{i}": float(p) for i, p in enumerate(pp)} | {"n_snps": int(len(m))}

def simulate_shared_signal(n=200, causal=60, shared=True, seed=1):
    """Synthetic demo/test data: one causal variant (shared or distinct)."""
    rng = np.random.default_rng(seed)
    snp = [f"rs{i}" for i in range(n)]
    def trait(c):
        se = np.full(n, 0.05)
        b = rng.normal(0, 0.05, n)
        b[c] += 0.4
        return pd.DataFrame({"snp": snp, "beta": b, "se": se})
    return trait(causal), trait(causal if shared else (causal + 40) % n)
