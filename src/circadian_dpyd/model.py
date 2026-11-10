"""Three-arm cohort state-transition model.

Arms:
A: empiric dosing
B: static CPIC genotyping
C: static genotyping plus circadian-adjusted phenotyping

First and recurrent Grade 3-4 toxicity events are tracked separately.
Recurrent events are excluded from primary estimates unless explicitly
parameterized. They are not assumed to be absent clinically.
"""
from __future__ import annotations

import numpy as np


FULL, RED, TOX_FIRST, POST, TOX_REC, DONE, DEAD = range(7)
ARMS = ("A", "B", "C")


def _cyc(p_cum, cycles):
    """Convert cumulative course risk to a constant per-cycle probability."""
    p = np.asarray(p_cum, dtype=float)
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p >= 1)):
        raise ValueError("Cumulative event risks must be finite and in [0, 1)")
    if cycles < 1:
        raise ValueError("Treatment cycles must be positive")
    return 1.0 - (1.0 - p) ** (1.0 / cycles)


def _trace(p, n_cycles, horizon, discount, p_first_full, p_first_red,
           start_reduced, p_recurrent, recurrence_estimable):
    """Trace one subgroup over a vector of parameter draws."""
    n = len(p_first_full)
    s = np.zeros((n, 7), dtype=float)
    s[:, RED if start_reduced else FULL] = 1.0

    pf = _cyc(p_first_full, n_cycles)
    pr = _cyc(p_first_red, n_cycles)

    p_recurrent = np.broadcast_to(
        np.asarray(p_recurrent, dtype=float), (n,)
    )
    if np.any(~np.isfinite(p_recurrent)) or np.any(
        (p_recurrent < 0) | (p_recurrent >= 1)
    ):
        raise ValueError("Recurrent event probabilities must be in [0, 1)")

    p_death = p["p_death_tox"]
    bg = 1.0 - (1.0 - p["bg_mortality_yr"]) ** (1.0 / 12.0)

    utility = np.stack([
        p["u_on_treatment"],
        p["u_on_treatment"],
        p["u_tox"],
        p["u_on_treatment"],
        p["u_tox"],
        p["u_done"],
        np.zeros(n),
    ], axis=1)

    qaly = np.zeros(n)
    first_events = np.zeros(n)
    recurrent_events = np.zeros(n)
    first_cost = np.zeros(n)
    recurrent_cost = np.zeros(n)

    for t in range(horizon):
        d = (1.0 + discount) ** (-t / 12.0)
        nxt = np.zeros_like(s)

        if t < n_cycles:
            # Patients without a prior event remain at their assigned dose
            # unless they experience their first Grade 3-4 toxicity event.
            new_first = pf * s[:, FULL] + pr * s[:, RED]
            nxt[:, TOX_FIRST] += new_first
            nxt[:, FULL] += (1.0 - pf) * s[:, FULL]
            nxt[:, RED] += (1.0 - pr) * s[:, RED]

            # Patients with a prior event are in the post-event reduced-dose
            # state and may experience a recurrent event during treatment.
            new_recurrent = p_recurrent * s[:, POST]
            nxt[:, TOX_REC] += new_recurrent
            nxt[:, POST] += (1.0 - p_recurrent) * s[:, POST]

            # Resolve toxicity episodes carried over from the previous cycle.
            prior_tox = s[:, TOX_FIRST] + s[:, TOX_REC]
            nxt[:, POST] += (1.0 - p_death) * prior_tox
            nxt[:, DEAD] += p_death * prior_tox

            first_events += new_first
            recurrent_events += new_recurrent
            first_cost += d * p["cost_tox_event"] * new_first

            if recurrence_estimable:
                recurrent_cost += (
                    d * p["cost_recurrent_event"] * new_recurrent
                )
        else:
            # After active treatment, survivors transition to completed.
            prior_tox = s[:, TOX_FIRST] + s[:, TOX_REC]
            nxt[:, DONE] += (
                s[:, FULL] + s[:, RED] + s[:, POST]
                + (1.0 - p_death) * prior_tox
                + s[:, DONE]
            )
            nxt[:, DEAD] += s[:, DEAD] + p_death * prior_tox

        # Apply background mortality to all living states after transitions.
        dying = bg[:, None] * nxt[:, :DEAD]
        nxt[:, :DEAD] -= dying
        nxt[:, DEAD] += dying.sum(axis=1)

        # Trapezoidal state-utility integration (half-cycle correction).
        qaly += (
            d * ((s + nxt) / 2.0 * utility).sum(axis=1) / 12.0
        )
        s = nxt

    if not recurrence_estimable:
        recurrent_events[:] = np.nan

    return qaly, first_cost + recurrent_cost, first_events, recurrent_events


def evaluate(p: dict, prev: dict, cfg: dict,
             recurrence_override=None) -> dict:
    """Evaluate all three arms for a parameter set or PSA sample.

    recurrence_override is a per-cycle recurrent event probability. Supplying
    it activates the exploratory recurrence pathway. The primary analysis
    leaves it unset and reports recurrence as not estimated.
    """
    meta = cfg["meta"]
    horizon = int(meta["horizon_months"])
    treatment_cycles = int(meta["treatment_cycles"])
    discount = float(meta["discount"])

    recurrence_estimable = recurrence_override is not None
    if recurrence_estimable:
        recurrence = recurrence_override
    else:
        recurrence = cfg.get("recurrence", {}).get("p_recurrent_cycle")

    if recurrence is None:
        recurrence_estimable = False
        recurrence = 0.0

    p = {
        k: np.asarray(v, dtype=float)
        for k, v in p.items()
        if v is not None
    }

    n = len(np.atleast_1d(p["p_tox_cum"]))
    p = {
        k: np.broadcast_to(v, (n,)).copy()
        for k, v in p.items()
    }

    if recurrence_estimable:
        p_recurrent = np.broadcast_to(
            np.asarray(recurrence, dtype=float), (n,)
        )
        if "cost_recurrent_event" not in p:
            raise ValueError(
                "Exploratory recurrence requires cost_recurrent_event"
            )
    else:
        p_recurrent = np.zeros(n)

    reclass_conditional = p["reclass_rate"]
    if np.any((reclass_conditional < 0) | (reclass_conditional > 1)):
        raise ValueError(
            "Conditional NM reclassification rate must be between 0 and 1"
        )
    reclass = prev["NM"] * reclass_conditional

    def reduced(risk):
        effect = p["dose_red_effect"]
        if np.any((effect < 0) | (effect > 1)):
            raise ValueError("dose_red_effect must be between 0 and 1")
        return p["p_tox_cum"] + (risk - p["p_tox_cum"]) * (1.0 - effect)

    groups = {
        "NM_rest": (
            prev["NM"] - reclass,
            p["p_tox_cum"],
            None,
        ),
        "circadian_reclassified": (
            reclass,
            p["p_tox_cum"] * p["rr_im"],
            "C",
        ),
        "IM": (
            np.full(n, prev["IM"]),
            p["p_tox_cum"] * p["rr_im"],
            "BC",
        ),
        "PM": (
            np.full(n, prev["PM"]),
            p["p_tox_cum"] * p["rr_pm"],
            "BC",
        ),
    }

    out = {
        arm: {
            "cost": np.zeros(n),
            "qaly": np.zeros(n),
            "first_events": np.zeros(n),
            "recurrent_events": np.zeros(n),
        }
        for arm in ARMS
    }

    for weight, first_risk, reduce_in in groups.values():
        first_risk = np.minimum(first_risk, 0.999999)
        reduced_risk = np.minimum(reduced(first_risk), 0.999999)

        for arm in ARMS:
            start_reduced = reduce_in is not None and arm in reduce_in
            q, c, first, recurrent = _trace(
                p,
                treatment_cycles,
                horizon,
                discount,
                first_risk,
                reduced_risk,
                start_reduced,
                p_recurrent,
                recurrence_estimable,
            )

            out[arm]["qaly"] += weight * q
            out[arm]["cost"] += weight * c
            out[arm]["first_events"] += weight * first

            if recurrence_estimable:
                out[arm]["recurrent_events"] += weight * recurrent
            else:
                out[arm]["recurrent_events"][:] = np.nan

    out["B"]["cost"] += p["cost_genotyping"]
    out["C"]["cost"] += p["cost_genotyping"] + p["cost_decision_support"]

    return out


def incremental(out: dict, wtp: float) -> dict:
    """Incremental Arm C versus Arm B outcomes, per patient."""
    d_cost = out["C"]["cost"] - out["B"]["cost"]
    d_qaly = out["C"]["qaly"] - out["B"]["qaly"]
    first_prevented = (
        out["B"]["first_events"] - out["C"]["first_events"]
    )
    recurrent_prevented = (
        out["B"]["recurrent_events"] - out["C"]["recurrent_events"]
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        icer = np.where(d_qaly != 0, d_cost / d_qaly, np.nan)

    return {
        "d_cost": d_cost,
        "d_qaly": d_qaly,
        "first_events_prevented": first_prevented,
        "recurrent_events_prevented": recurrent_prevented,
        "icer": icer,
        "nmb": wtp * d_qaly - d_cost,
        "dominant": (d_cost <= 0) & (d_qaly >= 0),
    }
