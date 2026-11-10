"""Load parameters and draw probabilistic samples from moment-specified distributions."""
from __future__ import annotations
import numpy as np, yaml
from scipy import stats

def load(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)

def _dist(spec: dict, base: float | None):
    t = spec["type"]
    if t == "fixed":
        return None
    if t == "beta":
        m, s = spec["mean"], spec["sd"]
        k = m * (1 - m) / s**2 - 1
        return stats.beta(m * k, (1 - m) * k)
    if t == "gamma":
        m, s = spec["mean"], spec["sd"]
        return stats.gamma(a=(m / s) ** 2, scale=s**2 / m)
    if t == "lognormal":
        m = spec.get("mean", base)
        s = spec.get("sd", spec.get("rel_sd", 0) * m)
        sig2 = np.log(1 + (s / m) ** 2)
        return stats.lognorm(s=np.sqrt(sig2), scale=m / np.sqrt(1 + (s / m) ** 2))
    if t == "uniform":
        return stats.uniform(spec["low"], spec["high"] - spec["low"])
    raise ValueError(t)

def base_values(cfg: dict) -> dict:
    return {k: v["value"] for k, v in cfg["parameters"].items()}

def sample(cfg: dict, base: dict, n: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    out = {}
    for k, v in cfg["parameters"].items():
        d = _dist(v["dist"], base[k])
        out[k] = np.full(n, base[k], float) if d is None else d.rvs(size=n, random_state=rng)
    return out

def interval(cfg: dict, base: dict, name: str, q=(0.025, 0.975)):
    """Low/high values for one-way sensitivity (distribution quantiles)."""
    d = _dist(cfg["parameters"][name]["dist"], base[name])
    return (base[name], base[name]) if d is None else tuple(d.ppf(q))
