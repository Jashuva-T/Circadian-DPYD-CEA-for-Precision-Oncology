"""Scenario-specific base-case, PSA, CEAC and sensitivity analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import model, params as P


def base_case(base, prev, cfg):
    out = model.evaluate(base, prev, cfg)
    inc = model.incremental(out, cfg["meta"]["wtp"])

    rows = []
    for arm in model.ARMS:
        rows.append({
            "arm": arm,
            "cost_eur_per_patient": float(out[arm]["cost"][0]),
            "qaly_per_patient": float(out[arm]["qaly"][0]),
            "first_toxicity_events_per_100": float(
                out[arm]["first_events"][0] * 100
            ),
        })

    summary = {
        key: float(np.atleast_1d(value)[0])
        for key, value in inc.items()
        if key != "recurrent_events_prevented"
    }
    summary["first_toxicity_events_prevented_per_100"] = (
        summary.pop("first_events_prevented") * 100
    )
    summary["incremental_cost_eur_per_100"] = summary["d_cost"] * 100
    summary["incremental_qalys_per_100"] = summary["d_qaly"] * 100
    summary["recurrent_events_estimable"] = False
    return pd.DataFrame(rows), summary


def psa(cfg, base, prev, wtp, n, seed, reclass_rate):
    draws = P.sample(cfg, base, n, seed)
    draws["reclass_rate"] = np.full(n, reclass_rate, dtype=float)

    inc = model.incremental(model.evaluate(draws, prev, cfg), wtp)
    return pd.DataFrame({
        "d_cost": inc["d_cost"],
        "d_qaly": inc["d_qaly"],
        "icer": inc["icer"],
        "nmb": inc["nmb"],
        "dominant": inc["dominant"],
        "first_events_prevented_per_100":
            inc["first_events_prevented"] * 100,
        "reclassifications_per_100_nm": reclass_rate * 100,
    })


def ceac(df, wtps):
    return pd.DataFrame({
        "wtp": wtps,
        "prob_ce": [
            float(np.mean(value * df["d_qaly"] - df["d_cost"] > 0))
            for value in wtps
        ],
    })


def tornado(cfg, base, prev, wtp):
    reference = model.incremental(
        model.evaluate(base, prev, cfg), wtp
    )["nmb"][0]

    rows = []
    for name, spec in cfg["parameters"].items():
        if spec["dist"]["type"] == "fixed":
            continue

        low, high = P.interval(cfg, base, name)
        values = []
        for value in (low, high):
            candidate = dict(base)
            candidate[name] = value
            values.append(
                model.incremental(
                    model.evaluate(candidate, prev, cfg), wtp
                )["nmb"][0]
            )

        rows.append({
            "parameter": name,
            "low_value": float(low),
            "high_value": float(high),
            "nmb_low": float(values[0]),
            "nmb_high": float(values[1]),
            "swing": float(abs(values[1] - values[0])),
        })

    return reference, pd.DataFrame(rows).sort_values(
        "swing", ascending=False
    )


def threshold(cfg, base, prev, wtp, grid=None):
    # Grid is conditional reclassification per 100 baseline NM patients.
    grid = np.arange(0.0, 50.01, 0.5) if grid is None else grid
    rows = []

    for value in grid:
        candidate = dict(base)
        candidate["reclass_rate"] = float(value) / 100
        inc = model.incremental(
            model.evaluate(candidate, prev, cfg), wtp
        )
        rows.append({
            "reclassifications_per_100_nm": float(value),
            "nmb_per_patient": float(inc["nmb"][0]),
            "icer": float(inc["icer"][0]),
            "incremental_cost_per_patient": float(inc["d_cost"][0]),
            "incremental_qaly_per_patient": float(inc["d_qaly"][0]),
        })

    frame = pd.DataFrame(rows)
    values = frame["nmb_per_patient"].to_numpy()
    crossings = np.where((values[:-1] <= 0) & (values[1:] > 0))[0]
    break_even = None

    if len(crossings):
        i = int(crossings[0])
        x0, x1 = grid[i], grid[i + 1]
        y0, y1 = values[i], values[i + 1]
        break_even = float(x0 - y0 * (x1 - x0) / (y1 - y0))

    return frame, break_even
