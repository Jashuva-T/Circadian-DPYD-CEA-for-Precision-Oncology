"""Publication-style figures (navy/blue/amber palette; no red/green; colour-blind safe)."""
from __future__ import annotations
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import circadian

NAVY, BLUE, AMBER, GREY = "#0B2545", "#1F6FA8", "#F2A900", "#7A8794"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titleweight": "bold", "figure.dpi": 150})
NOTE = "In silico model output; author-selected assumptions and parameter uncertainty are documented in config/params.yaml."

def _save(fig, out, name):
    fig.text(0.01, 0.005, NOTE, fontsize=7, color=GREY, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/figures/{name}.{ext}", dpi=300 if ext == "png" else None)
    plt.close(fig)

def fig_oscillation(cfg, out):
    f, zp = cfg["circadian"]["fold_change"], cfg["circadian"]["zt_peak"]
    zt = np.arange(0, 24.01, 0.25); y = circadian.relative_activity(zt, f, zp)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(zt, y, color=BLUE, lw=2.5)
    ax.axhline(1, color=GREY, ls=":", lw=1)
    ax.scatter([zp % 24, (zp + 12) % 24], [y.max(), y.min()], color=AMBER, zorder=3, s=50)
    ax.set(xlabel="Zeitgeber time (ZT)", ylabel="Relative DPD activity (24 h mean = 1)", xlim=(0, 24),
           title=f"Sinusoidal DPD activity, {f:g}-fold peak/trough")
    ax.set_xticks(range(0, 25, 4)); _save(fig, out, "fig1_circadian_oscillation")

def fig_reclass_sensitivity(cfg, out):
    f, zp = cfg["circadian"]["fold_change"], cfg["circadian"]["zt_peak"]
    taus = np.arange(0.8, 2.01, 0.05)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for (w, cv), c, ls in [(((8, 18), 0.0), BLUE, "-"), (((8, 18), 0.25), NAVY, "--"), (((0, 24), 0.0), AMBER, "-")]:
        r = [100 * circadian.nm_reclassification_rate(t, w, f, zp, cv) for t in taus]
        ax.plot(taus, r, color=c, ls=ls, lw=2, label=f"dosing ZT{w[0]}-{w[1]}, inter-individual CV {cv:g}")
    lo, hi = cfg["reclassification"]["per_100_range"]
    ax.axhspan(lo, hi, color=GREY, alpha=0.25, label="abstract-reported 8-12 per 100")
    ax.set(xlabel="NM/IM score threshold (activity score)", ylabel="NM patients reclassified per 100",
           title="Reclassification is assumption-dependent (illustrative)")
    ax.legend(fontsize=8, frameon=False); _save(fig, out, "fig2_reclassification_sensitivity")

def fig_ce_plane(df, wtp, out):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(df.d_qaly * 100, df.d_cost * 100, s=4, alpha=0.25, color=BLUE, rasterized=True)
    x = np.array([df.d_qaly.min(), df.d_qaly.max()]) * 100
    ax.plot(x, wtp * x, color=AMBER, lw=2, label=f"WTP EUR {wtp:,.0f}/QALY")
    ax.axhline(0, color=GREY, lw=0.8); ax.axvline(0, color=GREY, lw=0.8)
    ax.set(xlabel="Incremental QALYs per 100 patients", ylabel="Incremental cost per 100 patients (EUR)",
           title="Cost-effectiveness plane: Arm C vs Arm B")
    ax.legend(frameon=False); _save(fig, out, "fig3_ce_plane")

def fig_ceac(c, wtp, out):
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(c.wtp, c.prob_ce, color=BLUE, lw=2.5); ax.axvline(wtp, color=AMBER, ls="--")
    ax.set(xlabel="Willingness to pay (EUR/QALY)", ylabel="Probability cost-effective", ylim=(0, 1.02),
           title="Acceptability curve, Arm C vs B")
    _save(fig, out, "fig4_ceac")

def fig_tornado(ref, tor, out, top=12):
    t = tor.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for i, r in enumerate(t.itertuples()):
        ax.barh(i, r.nmb_low - ref, left=ref, color=BLUE, height=0.6)
        ax.barh(i, r.nmb_high - ref, left=ref, color=AMBER, height=0.6)
    ax.axvline(ref, color=NAVY, lw=1); ax.axvline(0, color=GREY, ls=":", lw=1)
    ax.set_yticks(range(len(t))); ax.set_yticklabels(t.parameter)
    ax.set(xlabel="Net monetary benefit per patient (EUR) at WTP",
           title="One-way sensitivity (blue = low, amber = high)")
    _save(fig, out, "fig5_tornado")

def fig_threshold(df, be, out):
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(df.reclassifications_per_100_nm, df.nmb_per_patient * 100, color=BLUE, lw=2.5); ax.axhline(0, color=GREY, lw=1)
    if be is not None:
        ax.axvline(be, color=AMBER, ls="--"); ax.text(be, ax.get_ylim()[1] * 0.9, f" break-even {be:.1f}/100", color=NAVY)
    ax.set(xlabel="NM-to-IM reclassifications per 100 patients", ylabel="Net monetary benefit per 100 (EUR)",
           title="Threshold analysis (effect per patient fixed)")
    _save(fig, out, "fig6_threshold")

def fig_psa_bars(summary, out):
    fig, ax = plt.subplots(figsize=(5, 3.8))
    v = [100 * summary["prob_ce"], 100 * summary["prob_dominant"]]
    b = ax.bar(["Cost-effective", "Dominant"], v, color=[BLUE, AMBER])
    ax.bar_label(b, fmt="%.1f%%"); ax.set(ylim=(0, 105), ylabel="% of PSA iterations", title="PSA outcomes (WTP EUR 30,000/QALY)")
    _save(fig, out, "fig7_psa_outcomes")


def fig_scenario_reclassification(results, out):
    """Compare mechanistic reclassification across the two dosing windows."""
    names = list(results)
    labels = [name.replace("_", " ").title() for name in names]
    values = [100 * results[name]["reclassification_rate"] for name in names]

    fig, ax = plt.subplots(figsize=(6.5, 4))
    bars = ax.bar(labels, values, color=[BLUE, AMBER][:len(values)])
    ax.bar_label(bars, fmt="%.2f")
    ax.set(
        ylabel="Reclassified baseline NM patients per 100",
        title="Mechanistic reclassification by dosing-time scenario",
        ylim=(0, max(values + [1]) * 1.25),
    )
    ax.text(
        0.01, -0.20,
        "Calculated from the specified sinusoidal activity model; "
        "not calibrated to abstract-reported reclassification.",
        transform=ax.transAxes, fontsize=8, color=GREY,
    )
    _save(fig, out, "fig2_scenario_reclassification")
