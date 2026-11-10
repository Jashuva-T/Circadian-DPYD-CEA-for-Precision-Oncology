"""Diplotype-to-activity-score distribution under Hardy-Weinberg equilibrium.

This is a simplified allele-frequency model. It assumes independent alleles
under HWE and does not model linkage disequilibrium or phased haplotypes.
"""
from __future__ import annotations

import itertools


def activity_score_distribution(cfg: dict) -> dict[float, float]:
    alleles = cfg["genetics"]["alleles"]

    frequencies = {
        name: float(spec["freq"])
        for name, spec in alleles.items()
    }
    activities = {
        name: float(spec["activity"])
        for name, spec in alleles.items()
    }

    if any(freq < 0 or freq > 1 for freq in frequencies.values()):
        raise ValueError("Allele frequencies must lie between 0 and 1")

    total_variant_frequency = sum(frequencies.values())
    if total_variant_frequency > 1:
        raise ValueError("Sum of variant allele frequencies cannot exceed 1")

    if any(activity < 0 for activity in activities.values()):
        raise ValueError("Allele activity values must be non-negative")

    frequencies["ref"] = 1.0 - total_variant_frequency
    activities["ref"] = 1.0

    distribution: dict[float, float] = {}

    for allele_a, allele_b in itertools.product(frequencies, repeat=2):
        score = activities[allele_a] + activities[allele_b]
        probability = frequencies[allele_a] * frequencies[allele_b]
        distribution[score] = distribution.get(score, 0.0) + probability

    if abs(sum(distribution.values()) - 1.0) > 1e-10:
        raise ArithmeticError("Diplotype probabilities do not sum to one")

    return dict(sorted(distribution.items()))


def phenotype_prevalence(cfg: dict) -> dict:
    distribution = activity_score_distribution(cfg)
    nm_threshold = float(cfg["genetics"]["nm_threshold"])
    im_threshold = float(cfg["genetics"].get("im_threshold", 1.0))

    if im_threshold < 0 or nm_threshold <= im_threshold:
        raise ValueError("Require 0 <= im_threshold < nm_threshold")

    prevalence = {"NM": 0.0, "IM": 0.0, "PM": 0.0}

    for score, probability in distribution.items():
        if score >= nm_threshold:
            group = "NM"
        elif score >= im_threshold:
            group = "IM"
        else:
            group = "PM"
        prevalence[group] += probability

    if abs(sum(prevalence.values()) - 1.0) > 1e-10:
        raise ArithmeticError("Phenotype prevalences do not sum to one")

    return prevalence
