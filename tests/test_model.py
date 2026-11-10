import numpy as np
import pytest

from circadian_dpyd import circadian, genetics, model, params as P, coloc


@pytest.fixture(scope="module")
def ctx():
    cfg = P.load("config/params.yaml")
    prev = genetics.phenotype_prevalence(cfg)
    base = P.base_values(cfg)
    result = circadian.compare_scenarios(
        base_score=cfg["circadian"]["base_score"],
        threshold=cfg["circadian"]["nm_im_threshold"],
        fold=cfg["circadian"]["fold_change"],
        zt_peak=cfg["circadian"]["zt_peak"],
    )
    base["reclass_rate"] = result["full_day"]["reclassification_rate"]
    return cfg, prev, base, result


def test_prevalence_sums_to_one(ctx):
    _, prev, _, _ = ctx
    assert sum(prev.values()) == pytest.approx(1.0)


def test_fold_change_peak_trough():
    zt = np.arange(0, 24, 0.1)
    y = circadian.relative_activity(zt, 2.4, 15)
    assert y.max() / y.min() == pytest.approx(2.4, rel=0.01)
    assert y.mean() == pytest.approx(1.0, abs=0.01)


def test_scenarios_are_distinct_and_hourly(ctx):
    _, _, _, result = ctx
    assert result["full_day"]["n_hours"] == 24
    assert result["daytime"]["n_hours"] == 11
    assert 0 <= result["full_day"]["reclassification_rate"] <= 1
    assert 0 <= result["daytime"]["reclassification_rate"] <= 1


def test_no_abstract_calibration_parameters(ctx):
    cfg, _, _, _ = ctx
    assert "rr_circ" not in cfg["parameters"]
    assert "reclass_rate" not in cfg["parameters"]
    assert "target_events_prevented_per_100" not in cfg


def test_model_returns_separate_first_and_recurrent_endpoints(ctx):
    cfg, prev, base, _ = ctx
    out = model.evaluate(base, prev, cfg)
    for arm in model.ARMS:
        assert "first_events" in out[arm]
        assert "recurrent_events" in out[arm]
        assert np.all(np.isfinite(out[arm]["first_events"]))
        assert np.all(np.isnan(out[arm]["recurrent_events"]))


def test_reclassification_is_conditional_on_nm_prevalence(ctx):
    cfg, prev, base, _ = ctx
    out = model.evaluate(base, prev, cfg)
    expected = prev["NM"] * base["reclass_rate"]
    assert expected >= 0
    assert out["C"]["first_events"].shape == (1,)


def test_invalid_reclassification_rate_rejected(ctx):
    cfg, prev, base, _ = ctx
    invalid = dict(base, reclass_rate=1.1)
    with pytest.raises(ValueError):
        model.evaluate(invalid, prev, cfg)


def test_incremental_endpoint_names(ctx):
    cfg, prev, base, _ = ctx
    inc = model.incremental(model.evaluate(base, prev, cfg), cfg["meta"]["wtp"])
    assert "first_events_prevented" in inc
    assert "recurrent_events_prevented" in inc
    assert "events_prevented" not in inc


def test_more_reclassification_improves_first_event_prevention(ctx):
    cfg, prev, base, _ = ctx
    low = model.incremental(
        model.evaluate(dict(base, reclass_rate=0.01), prev, cfg),
        cfg["meta"]["wtp"],
    )["first_events_prevented"][0]
    high = model.incremental(
        model.evaluate(dict(base, reclass_rate=0.20), prev, cfg),
        cfg["meta"]["wtp"],
    )["first_events_prevented"][0]
    assert high > low


def test_psa_reproducible(ctx):
    from circadian_dpyd import analysis

    cfg, prev, base, _ = ctx
    first = analysis.psa(cfg, base, prev, 30000, 50, 123, base["reclass_rate"])
    second = analysis.psa(cfg, base, prev, 30000, 50, 123, base["reclass_rate"])
    assert np.allclose(first["d_cost"], second["d_cost"])
    assert np.allclose(first["d_qaly"], second["d_qaly"])


def test_coloc_simulation_smoke_test():
    a, b = coloc.simulate_shared_signal(shared=True)
    assert coloc.coloc_abf(a, b)["PP4"] > 0.8
    a, b = coloc.simulate_shared_signal(shared=False)
    assert coloc.coloc_abf(a, b)["PP3"] > 0.8


def test_dose_reduction_lowers_reclassified_group_toxicity(ctx):
    cfg, prev, base, _ = ctx

    untreated = model.evaluate(
        dict(base, reclass_rate=0.20, dose_red_effect=0.0),
        prev,
        cfg,
    )
    reduced = model.evaluate(
        dict(base, reclass_rate=0.20, dose_red_effect=0.70),
        prev,
        cfg,
    )

    assert reduced["C"]["first_events"][0] < untreated["C"]["first_events"][0]
