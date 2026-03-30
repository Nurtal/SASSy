#!/usr/bin/env python3
"""Generate figures for all 5 OISA simulation runs."""

import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

RESULTS_DIR = Path(__file__).parent
FIGS_DIR = RESULTS_DIR / "figures"
FIGS_DIR.mkdir(exist_ok=True)

PALETTE = {
    "HSC":    "#2166ac",
    "CLP":    "#4dac26",
    "DN1":    "#d01c8b",
    "DN_pool":"#b35806",
    "DP":     "#542788",
    "CD4SP":  "#e08214",
    "CD8SP":  "#7fbc41",
    "CD4":    "#e08214",
    "CD8":    "#7fbc41",
    "pos_sel":"#4dac26",
    "neg_del":"#d73027",
    "neglect":"#878787",
    "export": "#1a9850",
    "flux":   "#2166ac",
    "lag":    "#d01c8b",
}

ENTITY_LABELS = {
    "CL:0000037": "HSC",
    "CL:0000051": "CLP",
    "CL:0002420": "DN1/ETP",
    "CL:0000893": "DP thymocyte",
    "CL:0000624": "CD4+ SP",
    "CL:0000625": "CD8+ SP",
    "CL:0000898": "Naïve T cell (export)",
    "CL:0000624_pln": "Naïve CD4+ T",
    "CL:0000625_pln": "Naïve CD8+ T",
}


# ─── helpers ──────────────────────────────────────────────────────────────────

def load_oissl(run_dir: Path) -> list[dict]:
    oissl_dir = run_dir / "oissl"
    files = sorted(oissl_dir.glob("checkpoint_*.json"),
                   key=lambda f: int(f.stem.split("_")[1]))
    return [json.loads(f.read_text()) for f in files]


def load_issl(run_dir: Path, model_id: str) -> list[dict]:
    issl_dir = run_dir / "issl" / model_id
    files = sorted(issl_dir.glob("checkpoint_*.json"),
                   key=lambda f: int(f.stem.split("_")[1]))
    return [json.loads(f.read_text()) for f in files]


def s_to_day(s):
    return np.asarray(s) / 86400.0


def get_entity_series(issl_records, entity_id):
    """Return (times_days, counts, ci_lo, ci_hi) for a given entity_id."""
    times, counts, ci_lo, ci_hi = [], [], [], []
    for rec in issl_records:
        t = rec["envelope"]["sim_time_s"]
        for ent in rec.get("continuous_state", []):
            if ent["entity_id"] == entity_id:
                times.append(t)
                counts.append(ent["count"])
                ci = ent.get("ci_95") or [ent["count"], ent["count"]]
                ci_lo.append(ci[0])
                ci_hi.append(ci[1])
                break
    return s_to_day(times), np.array(counts), np.array(ci_lo), np.array(ci_hi)


def get_export_flux(issl_records):
    """Return (times_days, fluxes) from export_signals."""
    times, fluxes = [], []
    for rec in issl_records:
        t = rec["envelope"]["sim_time_s"]
        for sig in rec.get("export_signals", []):
            times.append(t)
            fluxes.append(sig.get("flux", 0))
            break
    return s_to_day(times), np.array(fluxes)


def get_bio_export_flux(issl_records):
    """Return (times_days, bio_fluxes_per_day) reading biological_flux_per_day from export_signals.
    Falls back to flux * scaling_factor if biological_flux_per_day is absent."""
    times, fluxes = [], []
    for rec in issl_records:
        t = rec["envelope"]["sim_time_s"]
        for sig in rec.get("export_signals", []):
            times.append(t)
            bio = sig.get("biological_flux_per_day")
            if bio is None:
                sf = sig.get("scaling_factor", 1.0)
                bio = sig.get("flux", 0) * sf
            fluxes.append(bio)
            break
    return s_to_day(times), np.array(fluxes)


def get_discrete_events(issl_records):
    """Return dict {event_label: (times_days, counts)} from discrete_events."""
    data: dict[str, tuple] = {}
    for rec in issl_records:
        t = rec["envelope"]["sim_time_s"]
        for ev in rec.get("discrete_events", []):
            label = ev.get("label", ev.get("event_type", ev.get("event_id", "unknown")))
            if label not in data:
                data[label] = ([], [])
            data[label][0].append(t)
            data[label][1].append(ev.get("count", 0))
    return {k: (s_to_day(v[0]), np.array(v[1])) for k, v in data.items()}


def styled_ax(ax, title="", xlabel="Day", ylabel=""):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")


def fill_ci(ax, days, lo, hi, color, alpha=0.15):
    """Fill CI band, clamping negative lower bounds to 0."""
    lo_c = np.clip(lo, 0, None)
    ax.fill_between(days, lo_c, hi, color=color, alpha=alpha)


# ─── Fig 1 — BM1 baseline ─────────────────────────────────────────────────────

def fig_bm1_baseline(run_dir):
    bm = load_issl(run_dir, "bm_haematopoiesis")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    fig.suptitle("BM1 — Bone Marrow Haematopoiesis Baseline v7 (30 days, analytical SS init)", fontsize=12)

    entities = [
        ("CL:0000037", "HSC", PALETTE["HSC"]),
        ("CL:0000051", "CLP", PALETTE["CLP"]),
        ("CL:0002420", "DN1/ETP", PALETTE["DN1"]),
    ]

    for ax, (eid, label, color) in zip(axes, entities):
        days, counts, ci_lo, ci_hi = get_entity_series(bm, eid)
        ax.plot(days, counts, color=color, lw=2, label=label)
        fill_ci(ax, days, ci_lo, ci_hi, color)
        styled_ax(ax, title=label, ylabel="Cell count")
        ax.set_xlim(0, 30)

    # Inset: DN1 export flux with equilibrium reference line
    ax = axes[2]
    ax2 = ax.inset_axes([0.55, 0.55, 0.42, 0.38])
    t_flux, flux = get_export_flux(bm)
    ax2.plot(t_flux, flux, color=PALETTE["flux"], lw=1.5)
    ax2.axhline(59.5547, color="gray", ls="--", lw=0.8, alpha=0.7)   # v6 analytical SS
    ax2.set_title("Export flux", fontsize=7)
    ax2.set_xlabel("day", fontsize=6)
    ax2.set_ylabel("cells·day⁻¹", fontsize=6)
    ax2.tick_params(labelsize=6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.savefig(FIGS_DIR / "fig1_BM1_baseline.png", dpi=150)
    plt.close(fig)
    print("✓ fig1_BM1_baseline.png")


# ─── Fig 2 — THY1 baseline ────────────────────────────────────────────────────

def fig_thy1_baseline(run_dir):
    thy = load_issl(run_dir, "thymus_selection")

    fig = plt.figure(figsize=(13, 8), constrained_layout=True)
    fig.suptitle("THY1 — Thymic T-Cell Selection Baseline v3 (30 days, scaling_factor=300 000)", fontsize=12)
    gs = gridspec.GridSpec(2, 3, figure=fig)

    stage_info = [
        ("CL:0002420", "DN1/ETP",  PALETTE["DN1"],   gs[0, 0]),
        ("CL:0000893", "DP",       PALETTE["DP"],    gs[0, 1]),
        ("CL:0000624", "CD4+ SP",  PALETTE["CD4SP"], gs[0, 2]),
        ("CL:0000625", "CD8+ SP",  PALETTE["CD8SP"], gs[1, 0]),
    ]

    for eid, label, color, gspec in stage_info:
        ax = fig.add_subplot(gspec)
        days, counts, ci_lo, ci_hi = get_entity_series(thy, eid)
        ax.plot(days, counts / 1e6, color=color, lw=2)
        fill_ci(ax, days, ci_lo / 1e6, ci_hi / 1e6, color)
        styled_ax(ax, title=label, ylabel="Cells (×10⁶, scaled)")
        ax.set_xlim(0, 30)

    # Naïve T export — read from export_signals (not continuous_state)
    ax_exp = fig.add_subplot(gs[1, 1])
    t_exp, bio_flux = get_bio_export_flux(thy)
    ax_exp.fill_between(t_exp, bio_flux / 1e6, color=PALETTE["export"], alpha=0.4)
    ax_exp.plot(t_exp, bio_flux / 1e6, color=PALETTE["export"], lw=2)
    styled_ax(ax_exp, title="Naïve T export", ylabel="Cells·day⁻¹ (×10⁶, scaled)")
    ax_exp.set_xlim(0, 30)

    # Discrete events: selection counts
    ax_ev = fig.add_subplot(gs[1, 2])
    events = get_discrete_events(thy)
    ev_colors = {
        "positive_selection": PALETTE["pos_sel"],
        "negative_selection": PALETTE["neg_del"],
        "neglect_death":      PALETTE["neglect"],
    }
    plotted = False
    for ev_label, (t, c) in events.items():
        color = ev_colors.get(ev_label, "#555555")
        ax_ev.plot(t, c, lw=1.5, color=color, label=ev_label.replace("_", " ").title())
        plotted = True
    if plotted:
        ax_ev.legend(fontsize=7, frameon=False)
    styled_ax(ax_ev, title="Selection events / checkpoint", ylabel="Count")
    ax_ev.set_xlim(0, 30)

    fig.savefig(FIGS_DIR / "fig2_THY1_baseline.png", dpi=150)
    plt.close(fig)
    print("✓ fig2_THY1_baseline.png")


# ─── Fig 3 — COMP1 direct coupling ───────────────────────────────────────────

def fig_comp1_direct(run_dir):
    bm = load_issl(run_dir, "bm_haematopoiesis")
    thy = load_issl(run_dir, "thymus_selection")

    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    fig.suptitle("COMP1 — BM → Thymus Direct Coupling (no lag)", fontsize=12)

    # BM row
    for ax, (eid, label, color) in zip(axes[0], [
        ("CL:0000037", "HSC",     PALETTE["HSC"]),
        ("CL:0000051", "CLP",     PALETTE["CLP"]),
        ("CL:0002420", "DN1/ETP", PALETTE["DN1"]),
    ]):
        days, counts, ci_lo, ci_hi = get_entity_series(bm, eid)
        ax.plot(days, counts, color=color, lw=2)
        fill_ci(ax, days, ci_lo, ci_hi, color)
        styled_ax(ax, title=f"BM · {label}", ylabel="Cells")
        ax.set_xlim(0, 30)

    # Thymus row (counts are scaled by 300 000 — display in millions)
    for ax, (eid, label, color) in zip(axes[1], [
        ("CL:0000893", "DP",      PALETTE["DP"]),
        ("CL:0000624", "CD4+ SP", PALETTE["CD4SP"]),
        ("CL:0000625", "CD8+ SP", PALETTE["CD8SP"]),
    ]):
        days, counts, ci_lo, ci_hi = get_entity_series(thy, eid)
        ax.plot(days, counts / 1e6, color=color, lw=2)
        fill_ci(ax, days, ci_lo / 1e6, ci_hi / 1e6, color)
        styled_ax(ax, title=f"Thymus · {label}", ylabel="Cells (×10⁶, scaled)")
        ax.set_xlim(0, 30)

    # Annotate: DN1 export flux overlay on first thymus panel
    t_flux, flux = get_export_flux(bm)
    ax_twin = axes[1, 0].twinx()
    ax_twin.plot(t_flux, flux, color=PALETTE["flux"], lw=1, ls="--", alpha=0.7,
                 label="BM export flux")
    ax_twin.set_ylabel("Export flux (cells·day⁻¹)", fontsize=8, color=PALETTE["flux"])
    ax_twin.tick_params(labelcolor=PALETTE["flux"], labelsize=7)
    ax_twin.spines["top"].set_visible(False)

    fig.savefig(FIGS_DIR / "fig3_COMP1_direct.png", dpi=150)
    plt.close(fig)
    print("✓ fig3_COMP1_direct.png")


# ─── Fig 4 — COMP2 transfer model (lag) ──────────────────────────────────────

def fig_comp2_transfer(run_dir):
    bm  = load_issl(run_dir, "bm_haematopoiesis")
    thy = load_issl(run_dir, "thymus_selection")
    bt  = load_issl(run_dir, "blood_transit")

    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    fig.suptitle("COMP2 — BM → Blood Transit → Thymus (model-derived lag)", fontsize=12)

    # BM
    for ax, (eid, label, color) in zip(axes[0, :2], [
        ("CL:0000037", "HSC", PALETTE["HSC"]),
        ("CL:0002420", "DN1/ETP", PALETTE["DN1"]),
    ]):
        days, counts, ci_lo, ci_hi = get_entity_series(bm, eid)
        ax.plot(days, counts, color=color, lw=2)
        fill_ci(ax, days, ci_lo, ci_hi, color)
        styled_ax(ax, title=f"BM · {label}", ylabel="Cells")
        ax.set_xlim(0, 30)

    # Blood transit — cells delivered & lag
    ax_bt = axes[0, 2]
    t_bt, flux_bt = get_export_flux(bt)
    ax_bt.bar(t_bt, flux_bt, width=0.8, color=PALETTE["flux"], alpha=0.7,
              label="Cells delivered")
    styled_ax(ax_bt, title="Blood Transit · Cells delivered", ylabel="Cells")
    ax_bt.set_xlim(0, 30)

    # Extract lag_s from blood_transit ISSL
    lags_s = []
    t_lag  = []
    for rec in bt:
        t = rec["envelope"]["sim_time_s"]
        for sig in rec.get("export_signals", []):
            lag = sig.get("lag_s")
            if lag is not None:
                t_lag.append(t / 86400.0)
                lags_s.append(lag / 3600.0)  # hours
    ax_lag = ax_bt.twinx()
    ax_lag.plot(t_lag, lags_s, color=PALETTE["lag"], lw=1.5, ls="--",
                marker="o", ms=3, label="Transit lag (h)")
    ax_lag.set_ylabel("Transit lag (h)", fontsize=8, color=PALETTE["lag"])
    ax_lag.tick_params(labelcolor=PALETTE["lag"], labelsize=7)
    ax_lag.spines["top"].set_visible(False)

    # Thymus (counts scaled — display in millions)
    for ax, (eid, label, color) in zip(axes[1], [
        ("CL:0000893", "DP",      PALETTE["DP"]),
        ("CL:0000624", "CD4+ SP", PALETTE["CD4SP"]),
        ("CL:0000625", "CD8+ SP", PALETTE["CD8SP"]),
    ]):
        days, counts, ci_lo, ci_hi = get_entity_series(thy, eid)
        ax.plot(days, counts / 1e6, color=color, lw=2)
        fill_ci(ax, days, ci_lo / 1e6, ci_hi / 1e6, color)
        styled_ax(ax, title=f"Thymus · {label}", ylabel="Cells (×10⁶, scaled)")
        ax.set_xlim(0, 30)

    fig.savefig(FIGS_DIR / "fig4_COMP2_transfer.png", dpi=150)
    plt.close(fig)
    print("✓ fig4_COMP2_transfer.png")


# ─── Fig 5 — COMP3 full graph ─────────────────────────────────────────────────

def fig_comp3_full_graph(run_dir):
    bm  = load_issl(run_dir, "bm_haematopoiesis")
    thy = load_issl(run_dir, "thymus_selection")
    pln = load_issl(run_dir, "peripheral_ln")

    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    fig.suptitle("COMP3 — Full Immune Ontogeny Graph: BM → Blood → Thymus → Peripheral LN",
                 fontsize=12)
    gs = gridspec.GridSpec(3, 4, figure=fig)

    # Row 0: BM
    for col, (eid, label, color) in enumerate([
        ("CL:0000037", "HSC",     PALETTE["HSC"]),
        ("CL:0000051", "CLP",     PALETTE["CLP"]),
        ("CL:0002420", "DN1/ETP", PALETTE["DN1"]),
    ]):
        ax = fig.add_subplot(gs[0, col])
        days, counts, ci_lo, ci_hi = get_entity_series(bm, eid)
        ax.plot(days, counts, color=color, lw=2)
        fill_ci(ax, days, ci_lo, ci_hi, color)
        styled_ax(ax, title=f"BM · {label}", ylabel="Cells")
        ax.set_xlim(0, 30)

    # BM export flux
    ax_flux = fig.add_subplot(gs[0, 3])
    t_flux, flux = get_export_flux(bm)
    ax_flux.plot(t_flux, flux, color=PALETTE["flux"], lw=2)
    styled_ax(ax_flux, title="BM · DN1 export flux", ylabel="cells·day⁻¹")
    ax_flux.set_xlim(0, 30)

    # Row 1: Thymus stages (counts scaled — display in millions)
    for col, (eid, label, color) in enumerate([
        ("CL:0000893", "DP",       PALETTE["DP"]),
        ("CL:0000624", "CD4+ SP",  PALETTE["CD4SP"]),
        ("CL:0000625", "CD8+ SP",  PALETTE["CD8SP"]),
    ]):
        ax = fig.add_subplot(gs[1, col])
        days, counts, ci_lo, ci_hi = get_entity_series(thy, eid)
        ax.plot(days, counts / 1e6, color=color, lw=2)
        fill_ci(ax, days, ci_lo / 1e6, ci_hi / 1e6, color)
        styled_ax(ax, title=f"Thymus · {label}", ylabel="Cells (×10⁶, scaled)")
        ax.set_xlim(0, 30)

    # Thymus selection events
    ax_ev = fig.add_subplot(gs[1, 3])
    events = get_discrete_events(thy)
    ev_colors = {
        "positive_selection": PALETTE["pos_sel"],
        "negative_selection": PALETTE["neg_del"],
        "neglect_death":      PALETTE["neglect"],
    }
    for ev_label, (t, c) in events.items():
        color = ev_colors.get(ev_label, "#555555")
        ax_ev.plot(t, c, lw=1.5, color=color, label=ev_label.replace("_", " ").title())
    ax_ev.legend(fontsize=7, frameon=False)
    styled_ax(ax_ev, title="Thymus · Selection events", ylabel="Count/checkpoint")
    ax_ev.set_xlim(0, 30)

    # Row 2: Peripheral LN
    # CD4 and CD8 pools
    ax_pln = fig.add_subplot(gs[2, :2])
    pln_entities = [
        ("CL:0000624", "Naïve CD4+", PALETTE["CD4"]),
        ("CL:0000625", "Naïve CD8+", PALETTE["CD8"]),
    ]
    for eid, label, color in pln_entities:
        days, counts, ci_lo, ci_hi = get_entity_series(pln, eid)
        ax_pln.plot(days, counts, color=color, lw=2, label=label)
        fill_ci(ax_pln, days, ci_lo, ci_hi, color)
    ax_pln.legend(fontsize=8, frameon=False)
    styled_ax(ax_pln, title="Peripheral LN · Naïve T-cell pools", ylabel="Cells")
    ax_pln.set_xlim(0, 30)

    # CD4/CD8 ratio
    ax_ratio = fig.add_subplot(gs[2, 2])
    _, cd4, _, _ = get_entity_series(pln, "CL:0000624")
    days_cd8, cd8, _, _ = get_entity_series(pln, "CL:0000625")
    if len(cd4) and len(cd8):
        min_len = min(len(cd4), len(cd8))
        ratio = cd4[:min_len] / np.maximum(cd8[:min_len], 1)
        ax_ratio.plot(days_cd8[:min_len], ratio, color="#b2182b", lw=2)
        ax_ratio.axhline(2.0, ls="--", lw=1, color="gray", label="Expected ~2:1")
        ax_ratio.legend(fontsize=7, frameon=False)
    styled_ax(ax_ratio, title="PLN · CD4/CD8 ratio", ylabel="Ratio")
    ax_ratio.set_xlim(0, 30)
    ax_ratio.set_ylim(0, 3.0)

    # Naïve T export from thymus
    ax_texport = fig.add_subplot(gs[2, 3])
    t_te, flux_te = get_export_flux(thy)
    ax_texport.fill_between(t_te, flux_te, color=PALETTE["export"], alpha=0.6)
    ax_texport.plot(t_te, flux_te, color=PALETTE["export"], lw=1.5)
    styled_ax(ax_texport, title="Thymus → PLN export", ylabel="cells·checkpoint⁻¹")
    ax_texport.set_xlim(0, 30)

    fig.savefig(FIGS_DIR / "fig5_COMP3_full_graph.png", dpi=150)
    plt.close(fig)
    print("✓ fig5_COMP3_full_graph.png")


# ─── Fig 6 — comparative summary ─────────────────────────────────────────────

def fig_comparative(runs):
    """DN1 export flux comparison across COMP1, COMP2, COMP3."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    fig.suptitle("Comparative — DN1 Export Flux and Thymic Output Across Runs", fontsize=12)

    run_styles = {
        "COMP1_direct":   {"color": "#1b7837", "ls": "-",  "label": "COMP1 (direct)"},
        "COMP2_transfer": {"color": "#762a83", "ls": "--", "label": "COMP2 (transfer)"},
        "COMP3_full_graph": {"color": "#d6604d", "ls": "-.", "label": "COMP3 (full graph)"},
    }

    for run_id, run_dir in runs.items():
        style = run_styles.get(run_id)
        if style is None:
            continue
        try:
            bm  = load_issl(run_dir, "bm_haematopoiesis")
            thy = load_issl(run_dir, "thymus_selection")
        except Exception:
            continue

        t_bm, flux_bm = get_export_flux(bm)
        axes[0].plot(t_bm, flux_bm, **{k: v for k, v in style.items()}, lw=2)

        t_thy, flux_thy = get_export_flux(thy)
        axes[1].plot(t_thy, flux_thy, **{k: v for k, v in style.items()}, lw=2)

    styled_ax(axes[0], title="BM DN1 export flux", ylabel="cells·day⁻¹")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].set_xlim(0, 30)
    axes[0].set_ylim(bottom=0)

    styled_ax(axes[1], title="Thymic naïve T export flux", ylabel="cells·checkpoint⁻¹")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].set_xlim(0, 30)

    fig.savefig(FIGS_DIR / "fig6_comparative.png", dpi=150)
    plt.close(fig)
    print("✓ fig6_comparative.png")


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    run_dirs = {
        "BM1_baseline":    RESULTS_DIR / "BM1_baseline"  / "BM1_baseline",
        "THY1_baseline":   RESULTS_DIR / "THY1_baseline" / "THY1_baseline",
        "COMP1_direct":    RESULTS_DIR / "COMP1_direct"  / "COMP1_direct",
        "COMP2_transfer":  RESULTS_DIR / "COMP2_transfer"/ "COMP2_transfer",
        "COMP3_full_graph":RESULTS_DIR / "COMP3_full_graph"/"COMP3_full_graph",
    }

    print("Generating figures …")
    fig_bm1_baseline(run_dirs["BM1_baseline"])
    fig_thy1_baseline(run_dirs["THY1_baseline"])
    fig_comp1_direct(run_dirs["COMP1_direct"])
    fig_comp2_transfer(run_dirs["COMP2_transfer"])
    fig_comp3_full_graph(run_dirs["COMP3_full_graph"])
    fig_comparative(run_dirs)
    print(f"\nAll figures saved to {FIGS_DIR}/")


if __name__ == "__main__":
    main()
