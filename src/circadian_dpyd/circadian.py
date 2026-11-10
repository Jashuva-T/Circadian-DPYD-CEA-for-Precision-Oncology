"""Hourly circadian DPYD activity and scenario-specific reclassification.

The peak-to-trough expression/activity fold change is a model input, not
independent evidence of clinical DPD enzyme activity. The expression-to-
activity relationship and dosing-time distribution remain assumptions.
"""
from __future__ import annotations

import numpy as np


SCENARIOS = {
    "full_day": tuple(range(24)),
    "daytime": tuple(range(8, 19)),
}


def amplitude(fold: float) -> float:
    if not np.isfinite(fold) or fold < 1:
        raise ValueError("fold must be finite and >= 1")
    return (fold - 1.0) / (fold + 1.0)


def relative_activity(zt, fold: float = 2.4, zt_peak: float = 15.0):
    zt = np.asarray(zt, dtype=float)
    if not np.all(np.isfinite(zt)):
        raise ValueError("ZT values must be finite")
    return 1.0 + amplitude(fold) * np.cos(
        2.0 * np.pi * (zt - zt_peak) / 24.0
    )


def adjusted_score(base_score, zt, fold: float = 2.4, zt_peak: float = 15.0):
    return np.asarray(base_score, dtype=float) * relative_activity(
        zt, fold, zt_peak
    )


def hourly_profile(
    base_score: float = 2.0,
    threshold: float = 1.5,
    fold: float = 2.4,
    zt_peak: float = 15.0,
) -> list[dict]:
    """Return the 24-hour activity-score profile and threshold classification."""
    if base_score < 0 or threshold < 0:
        raise ValueError("Activity score and threshold must be non-negative")

    hours = np.arange(24)
    scores = adjusted_score(base_score, hours, fold, zt_peak)

    return [
        {
            "zt": int(hour),
            "adjusted_score": float(score),
            "reclassified": bool(score < threshold),
        }
        for hour, score in zip(hours, scores)
    ]


def scenario_reclassification(
    window,
    base_score: float = 2.0,
    threshold: float = 1.5,
    fold: float = 2.4,
    zt_peak: float = 15.0,
) -> dict:
    """Calculate reclassification under a discrete, uniform dosing-time window."""
    hours = np.asarray(tuple(window), dtype=int)

    if hours.size == 0:
        raise ValueError("Dosing-time window cannot be empty")
    if np.any((hours < 0) | (hours > 23)):
        raise ValueError("Dosing times must be integer ZT values from 0 to 23")
    if np.unique(hours).size != hours.size:
        raise ValueError("Dosing-time window cannot contain duplicate hours")

    scores = adjusted_score(base_score, hours, fold, zt_peak)
    classified = scores < threshold

    return {
        "n_hours": int(hours.size),
        "hours": hours.tolist(),
        "reclassified_hours": hours[classified].tolist(),
        "reclassification_rate": float(np.mean(classified)),
        "reclassifications_per_100": float(np.mean(classified) * 100),
        "mean_adjusted_score": float(np.mean(scores)),
        "minimum_adjusted_score": float(np.min(scores)),
        "maximum_adjusted_score": float(np.max(scores)),
    }


def compare_scenarios(
    base_score: float = 2.0,
    threshold: float = 1.5,
    fold: float = 2.4,
    zt_peak: float = 15.0,
) -> dict:
    """Return independent full-day and daytime scenario results."""
    return {
        name: scenario_reclassification(
            hours, base_score, threshold, fold, zt_peak
        )
        for name, hours in SCENARIOS.items()
    }


def nm_reclassification_rate(
    threshold=1.5,
    window=(8, 18),
    fold=2.4,
    zt_peak=15.0,
    cv=0.0,
    n=200_000,
    seed=7,
):
    """Compatibility function for the earlier continuous-time illustration.

    The primary analysis uses discrete hourly scenarios through
    compare_scenarios(), not this continuous-time approximation.
    """
    if cv < 0:
        raise ValueError("cv must be non-negative")
    if n < 1:
        raise ValueError("n must be positive")

    rng = np.random.default_rng(seed)
    zt = rng.uniform(float(window[0]), float(window[1]), n)
    spread = (
        np.exp(rng.normal(0, np.sqrt(np.log1p(cv**2)), n))
        if cv > 0 else 1.0
    )
    scores = adjusted_score(2.0, zt, fold, zt_peak) * spread
    return float(np.mean(scores < threshold))
