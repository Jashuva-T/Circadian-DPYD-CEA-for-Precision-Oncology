"""Reproducible two-scenario analysis pipeline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from . import analysis, circadian, figures, genetics, params as P


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/params.yaml")
    parser.add_argument("--out", default="results")
    args = parser.parse_args(argv)

    cfg = P.load(args.config)
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "figures").mkdir(parents=True, exist_ok=True)
    prev = genetics.phenotype_prevalence(cfg)
    base = P.base_values(cfg)

    circ_cfg = cfg["circadian"]
    scenario_results = {}
    reclassification = {}

    for name, hours in circ_cfg["scenarios"].items():
        result = circadian.scenario_reclassification(
            hours,
            base_score=circ_cfg["base_score"],
            threshold=circ_cfg["nm_im_threshold"],
            fold=circ_cfg["fold_change"],
            zt_peak=circ_cfg["zt_peak"],
        )
        reclassification[name] = result

        scenario_rate = float(result["reclassification_rate"])
        scenario_base = dict(base, reclass_rate=scenario_rate)
        scenario_dir = out_root / "scenarios" / name
        (scenario_dir / "tables").mkdir(parents=True, exist_ok=True)
        (scenario_dir / "figures").mkdir(parents=True, exist_ok=True)

        arms, base_summary = analysis.base_case(
            scenario_base, prev, cfg
        )
        arms.to_csv(
            scenario_dir / "tables" / "base_case_arms.csv",
            index=False,
        )

        draws = analysis.psa(
            cfg, base, prev, cfg["meta"]["wtp"],
            cfg["meta"]["n_psa"], cfg["meta"]["seed"],
            scenario_rate,
        )
        draws.to_csv(
            scenario_dir / "tables" / "psa_draws.csv.gz",
            index=False,
            compression="gzip",
            float_format="%.8g",
        )

        wtps = np.arange(0, 60001, 1000)
        acceptability = analysis.ceac(draws, wtps)
        acceptability.to_csv(
            scenario_dir / "tables" / "ceac.csv", index=False
        )

        # One-way sensitivity analysis for this scenario.
        reference, tornado = analysis.tornado(
            cfg, scenario_base, prev, cfg["meta"]["wtp"]
        )
        tornado.to_csv(
            scenario_dir / "tables" / "tornado.csv", index=False
        )

        threshold, break_even = analysis.threshold(
            cfg, scenario_base, prev, cfg["meta"]["wtp"]
        )
        threshold.to_csv(
            scenario_dir / "tables" / "threshold.csv", index=False
        )

        prob_ce = float(np.mean(draws["nmb"] > 0))
        prob_dominant = float(np.mean(draws["dominant"]))
        first_event_ci = [
            float(x) for x in np.nanpercentile(
                draws["first_events_prevented_per_100"], [2.5, 97.5]
            )
        ]

        scenario_summary = {
            "scenario": name,
            "dosing_hours_zt": list(hours),
            "reclassification": result,
            "base_case": base_summary,
            "probability_cost_effective": prob_ce,
            "probability_dominant": prob_dominant,
            "median_incremental_nmb_eur_per_patient": float(
                draws["nmb"].median()
            ),
            "median_incremental_cost_eur_per_patient": float(
                draws["d_cost"].median()
            ),
            "median_incremental_qaly_per_patient": float(
                draws["d_qaly"].median()
            ),
            "icer_interpretation": (
                "Not summarized: ICERs are not meaningful in dominant "
                "iterations, and are undefined when incremental QALYs are zero. "
                "See PSA draws for the iteration-level ICER values."
            ),
            "first_toxicity_events_prevented_per_100_95ci": first_event_ci,
            "break_even_reclassifications_per_100_nm": break_even,
            "recurrence": "Not estimated in primary analysis",
            "cost_per_first_event": (
                "Exploratory author-selected assumption; not independently "
                "sourced empirical unit cost"
            ),
        }
        scenario_results[name] = scenario_summary

        figures.fig_ce_plane(draws, cfg["meta"]["wtp"], str(scenario_dir))
        figures.fig_ceac(acceptability, cfg["meta"]["wtp"], str(scenario_dir))
        figures.fig_tornado(reference, tornado, str(scenario_dir))
        figures.fig_threshold(threshold, break_even, str(scenario_dir))

    figures.fig_oscillation(cfg, str(out_root))
    figures.fig_scenario_reclassification(reclassification, str(out_root))

    summary = {
        "model": "Three-arm cohort state-transition model",
        "arms": {
            "A": "Empiric dosing",
            "B": "Static CPIC DPYD genotyping",
            "C": "Static genotyping plus circadian-adjusted phenotyping",
        },
        "phenotype_prevalence": prev,
        "circadian_scenarios": reclassification,
        "scenario_results": scenario_results,
        "abstract_benchmarks_not_used_for_calibration": cfg[
            "abstract_benchmarks"
        ],
        "model_scope": {
            "horizon_months": cfg["meta"]["horizon_months"],
            "treatment_cycles": cfg["meta"]["treatment_cycles"],
            "discount_rate": cfg["meta"]["discount"],
            "perspective": "European payer scenario; cost inputs require validation",
            "recurrent_toxicity": "Not estimated in primary analysis",
            "clinical_validation": "Not performed by this computational model",
            "colocalization_pp4": (
                "Supplied abstract-reported value; not recomputed by this pipeline"
            ),
        },
        "interpretation": (
            "Outputs are conditional on author-selected assumptions. "
            "The two dosing-time scenarios are reported independently. "
            "Mechanistic reclassification is not forced to match abstract "
            "benchmarks. First toxicity events are reported separately from "
            "recurrent events, which are not estimated."
        ),
    }

    (out_root / "summary.json").write_text(
        json.dumps(_json_safe(summary), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(_json_safe({
        "phenotype_prevalence": prev,
        "scenario_results": {
            name: {
                "reclassification_per_100_nm":
                    item["reclassification"]["reclassifications_per_100"],
                "probability_cost_effective":
                    item["probability_cost_effective"],
                "break_even_reclassifications_per_100_nm":
                    item["break_even_reclassifications_per_100_nm"],
            }
            for name, item in scenario_results.items()
        },
    }), indent=2))


if __name__ == "__main__":
    main()
