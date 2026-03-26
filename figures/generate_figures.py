#!/usr/bin/env python3
"""
Generate figures for all 5 OISA simulation runs.

Usage:
    python figures/generate_figures.py

Outputs (saved to figures/):
    fig1_BM1_haematopoiesis.png   — BM baseline: HSC/CLP/DN1 dynamics + export flux
    fig2_THY1_selection.png       — Thymus baseline: stage populations + selection events
    fig3_COMP1_direct.png         — BM→Thymus direct coupling (zero lag)
    fig4_COMP2_transfer.png       — BM→Thymus via blood transit (biologically realistic lag)
    fig5_COMP3_full_graph.png     — Full pipeline: BM → transit → thymus → peripheral LN
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
FIGURES.mkdir(exist_ok=True)

# ── Style ─────────────────────────────────────────────────────────────────────
COLORS = {
    "HSC":   "#1f77b4",
    "MPP":   "#aec7e8",
    "LMPP":  "#6baed6",
    "CLP":   "#ff7f0e",
    "DN1":   "#2ca02c",
    "DN":    "#9467bd",
    "DP":    "#8c564b",
    "CD4SP": "#e377c2",
    "CD8SP": "#7f7f7f",
    "CD4":   "#d62728",
    "CD8":   "#17becf",
    "flux":  "#2ca02c",
    "pos_sel":  "#2ca02c",
    "neg_sel":  "#d62728",
    "neglect":  "#ff7f0e",
}

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "figure.dpi": 150,
})


# ── Loaders ───────────────────────────────────────────────────────────────────

def _sort_key(f: str) -> float:
    return float(Path(f).stem.split("_")[-1])


def load_bm(run_id: str = "latest") -> dict:
    """Load BM haematopoiesis checkpoints (v4: 5-compartment)."""
    files = sorted(
        glob.glob(str(RESULTS / run_id / "latest" / "issl" / "bm_haematopoiesis" / "checkpoint_*.json")),
        key=_sort_key,
    )
    days = []
    hsc, mpp, lmpp, clp, dn1 = [], [], [], [], []
    hsc_lo, hsc_hi = [], []
    mpp_lo, mpp_hi = [], []
    lmpp_lo, lmpp_hi = [], []
    clp_lo, clp_hi = [], []
    dn1_lo, dn1_hi = [], []
    flux, flux_lo, flux_hi = [], [], []

    for f in files:
        d = json.loads(Path(f).read_text())
        t = d["envelope"]["sim_time_s"] / 86400
        days.append(t)
        for ent in d["continuous_state"]:
            lbl = ent["label"].lower()
            c, lo, hi = ent["count"], ent["ci_95"][0], ent["ci_95"][1]
            if "haematopoietic stem" in lbl:
                hsc.append(c); hsc_lo.append(lo); hsc_hi.append(hi)
            elif "multipotent progenitor" in lbl and "lymphoid" not in lbl:
                mpp.append(c); mpp_lo.append(lo); mpp_hi.append(hi)
            elif "lymphoid-primed" in lbl:
                lmpp.append(c); lmpp_lo.append(lo); lmpp_hi.append(hi)
            elif "common lymphoid" in lbl:
                clp.append(c); clp_lo.append(lo); clp_hi.append(hi)
            elif "early t-lineage" in lbl or "dn1" in lbl:
                dn1.append(c); dn1_lo.append(lo); dn1_hi.append(hi)
        for sig in d.get("export_signals", []):
            flux.append(sig["flux"]); flux_lo.append(sig["ci_95"][0]); flux_hi.append(sig["ci_95"][1])
        if not d.get("export_signals"):
            flux.append(0); flux_lo.append(0); flux_hi.append(0)

    return dict(
        days=np.array(days),
        hsc=np.array(hsc), mpp=np.array(mpp), lmpp=np.array(lmpp),
        clp=np.array(clp), dn1=np.array(dn1),
        hsc_lo=np.array(hsc_lo), hsc_hi=np.array(hsc_hi),
        mpp_lo=np.array(mpp_lo), mpp_hi=np.array(mpp_hi),
        lmpp_lo=np.array(lmpp_lo), lmpp_hi=np.array(lmpp_hi),
        clp_lo=np.array(clp_lo), clp_hi=np.array(clp_hi),
        dn1_lo=np.array(dn1_lo), dn1_hi=np.array(dn1_hi),
        flux=np.array(flux), flux_lo=np.array(flux_lo), flux_hi=np.array(flux_hi),
    )


def load_thymus(run_prefix: str) -> dict:
    """Load thymus_selection checkpoints from a given results/*/latest path."""
    files = sorted(
        glob.glob(str(RESULTS / run_prefix / "latest" / "issl" / "thymus_selection" / "checkpoint_*.json")),
        key=_sort_key,
    )
    days = []
    dn, dp, cd4sp, cd8sp = [], [], [], []
    dn_lo, dn_hi = [], []
    dp_lo, dp_hi = [], []
    cd4sp_lo, cd4sp_hi = [], []
    cd8sp_lo, cd8sp_hi = [], []
    pos_sel, neg_sel, neglect = [], [], []
    export, export_lo, export_hi = [], [], []

    for f in files:
        d = json.loads(Path(f).read_text())
        t = d["envelope"]["sim_time_s"] / 86400
        days.append(t)

        for ent in d["continuous_state"]:
            lbl = ent["label"]
            if "Double-negative" in lbl:
                dn.append(ent["count"]); dn_lo.append(ent["ci_95"][0]); dn_hi.append(ent["ci_95"][1])
            elif "Double-positive" in lbl:
                dp.append(ent["count"]); dp_lo.append(ent["ci_95"][0]); dp_hi.append(ent["ci_95"][1])
            elif "CD4+" in lbl and "single" in lbl.lower():
                cd4sp.append(ent["count"]); cd4sp_lo.append(ent["ci_95"][0]); cd4sp_hi.append(ent["ci_95"][1])
            elif "CD8+" in lbl and "single" in lbl.lower():
                cd8sp.append(ent["count"]); cd8sp_lo.append(ent["ci_95"][0]); cd8sp_hi.append(ent["ci_95"][1])

        for ev in d.get("discrete_events", []):
            if ev["event_type"] == "positive_selection":
                pos_sel.append(ev["count"])
            elif ev["event_type"] == "negative_selection":
                neg_sel.append(ev["count"])
            elif ev["event_type"] == "neglect_death":
                neglect.append(ev["count"])

        sigs = d.get("export_signals", [])
        if sigs:
            export.append(sigs[0]["flux"])
            export_lo.append(sigs[0]["ci_95"][0])
            export_hi.append(sigs[0]["ci_95"][1])
        else:
            export.append(0); export_lo.append(0); export_hi.append(0)

    # Pad event lists if shorter (some checkpoints may lack events)
    n = len(days)
    for lst in [pos_sel, neg_sel, neglect]:
        while len(lst) < n:
            lst.append(0)

    return dict(
        days=np.array(days),
        dn=np.array(dn), dp=np.array(dp), cd4sp=np.array(cd4sp), cd8sp=np.array(cd8sp),
        dn_lo=np.array(dn_lo), dn_hi=np.array(dn_hi),
        dp_lo=np.array(dp_lo), dp_hi=np.array(dp_hi),
        cd4sp_lo=np.array(cd4sp_lo), cd4sp_hi=np.array(cd4sp_hi),
        cd8sp_lo=np.array(cd8sp_lo), cd8sp_hi=np.array(cd8sp_hi),
        pos_sel=np.array(pos_sel), neg_sel=np.array(neg_sel), neglect=np.array(neglect),
        export=np.array(export), export_lo=np.array(export_lo), export_hi=np.array(export_hi),
    )


def load_periph(run_prefix: str) -> dict:
    """Load peripheral_ln checkpoints."""
    files = sorted(
        glob.glob(str(RESULTS / run_prefix / "latest" / "issl" / "peripheral_ln" / "checkpoint_*.json")),
        key=_sort_key,
    )
    days, cd4, cd8 = [], [], []
    for f in files:
        d = json.loads(Path(f).read_text())
        days.append(d["envelope"]["sim_time_s"] / 86400)
        for ent in d["continuous_state"]:
            if "CD4" in ent["label"]:
                cd4.append(ent["count"])
            elif "CD8" in ent["label"]:
                cd8.append(ent["count"])
    return dict(days=np.array(days), cd4=np.array(cd4), cd8=np.array(cd8))


# ── Figure helpers ────────────────────────────────────────────────────────────

def _shade(ax, days, lo, hi, color, alpha=0.18):
    ax.fill_between(days, lo, hi, color=color, alpha=alpha, linewidth=0)


# ── Fig 1: BM1 ────────────────────────────────────────────────────────────────

def fig1_bm1():
    bm = load_bm("BM1")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1.5]})
    fig.suptitle("Fig 1 — BM1: Bone Marrow Haematopoiesis (ODE, 5-compartment, v4)", fontweight="bold")

    # Top: all 5 compartments (log scale)
    for key, label, color in [
        ("hsc",  "HSC",      COLORS["HSC"]),
        ("mpp",  "MPP",      COLORS["MPP"]),
        ("lmpp", "LMPP",     COLORS["LMPP"]),
        ("clp",  "CLP",      COLORS["CLP"]),
        ("dn1",  "DN1/ETP",  COLORS["DN1"]),
    ]:
        ax1.plot(bm["days"], bm[key], color=color, lw=1.8, label=label)
        _shade(ax1, bm["days"], bm[f"{key}_lo"], bm[f"{key}_hi"], color)

    ax1.set_ylabel("Cell count")
    ax1.set_yscale("log")
    ax1.legend(loc="upper right", ncol=2)
    ax1.set_title("Progenitor pool dynamics: HSC → MPP → LMPP → CLP → DN1 (log scale)", fontsize=9)
    ax1.grid(True, which="both", ls=":", alpha=0.4)

    # Annotate equilibrium values
    eq = {"HSC": 9075, "MPP": 1361, "LMPP": 595, "CLP": 340, "DN1": 181}
    for (key, lbl, col), (_, eq_val) in zip(
        [("hsc","HSC",COLORS["HSC"]),("mpp","MPP",COLORS["MPP"]),
         ("lmpp","LMPP",COLORS["LMPP"]),("clp","CLP",COLORS["CLP"]),
         ("dn1","ETP",COLORS["DN1"])],
        eq.items()
    ):
        ax1.axhline(eq_val, color=col, ls="--", lw=0.6, alpha=0.5)

    # Bottom: export flux
    ax2.plot(bm["days"], bm["flux"], color=COLORS["flux"], lw=1.8, label="DN1 export flux")
    _shade(ax2, bm["days"], bm["flux_lo"], bm["flux_hi"], COLORS["flux"])
    ax2.axhline(27.2, color="gray", ls="--", lw=0.8, alpha=0.6, label="Expected eq. (27.2 cells/day)")
    ax2.set_xlabel("Simulation time (days)")
    ax2.set_ylabel("cells · day⁻¹")
    ax2.set_title("Progenitor export flux to blood / thymus", fontsize=9)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.legend(loc="upper right")

    fig.tight_layout()
    out = FIGURES / "fig1_BM1_haematopoiesis.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Fig 2: THY1 ───────────────────────────────────────────────────────────────

def fig2_thy1():
    th = load_thymus("THY1")
    fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True,
                              gridspec_kw={"height_ratios": [2.5, 2, 1]})
    fig.suptitle("Fig 2 — THY1: Thymic T-cell Selection (ABM baseline)", fontweight="bold")

    ax1, ax2, ax3 = axes

    # Panel 1: stage populations
    for key, label, color in [
        ("dn",    "DN (all)",  COLORS["DN"]),
        ("dp",    "DP",        COLORS["DP"]),
        ("cd4sp", "CD4SP",     COLORS["CD4SP"]),
        ("cd8sp", "CD8SP",     COLORS["CD8SP"]),
    ]:
        ax1.plot(th["days"], th[key], color=color, lw=1.8, label=label)
        _shade(ax1, th["days"], th[f"{key}_lo"], th[f"{key}_hi"], color)
    ax1.set_ylabel("Mean cell count\n(12 realisations)")
    ax1.set_title("Thymocyte stage populations", fontsize=9)
    ax1.legend(loc="upper right", ncol=2)
    ax1.grid(True, ls=":", alpha=0.4)

    # Panel 2: selection events (stacked bar)
    width = np.diff(th["days"]).mean() * 0.8 if len(th["days"]) > 1 else 0.8
    ax2.bar(th["days"], th["neglect"],  width=width, label="Neglect death",    color=COLORS["neglect"],  alpha=0.8)
    ax2.bar(th["days"], th["pos_sel"], width=width, label="Positive selection",color=COLORS["pos_sel"],  alpha=0.8,
            bottom=th["neglect"])
    ax2.bar(th["days"], th["neg_sel"], width=width, label="Negative selection",color=COLORS["neg_sel"],  alpha=0.8,
            bottom=th["neglect"] + th["pos_sel"])
    ax2.set_ylabel("Events · checkpoint⁻¹")
    ax2.set_title("Selection events (mean across realisations)", fontsize=9)
    ax2.legend(loc="upper right", ncol=3)
    ax2.grid(True, axis="y", ls=":", alpha=0.4)

    # Panel 3: export flux
    ax3.plot(th["days"], th["export"], color=COLORS["flux"], lw=1.8)
    _shade(ax3, th["days"], th["export_lo"], th["export_hi"], COLORS["flux"])
    ax3.set_xlabel("Simulation time (days)")
    ax3.set_ylabel("cells · day⁻¹")
    ax3.set_title("Naïve T cell export", fontsize=9)
    ax3.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    out = FIGURES / "fig2_THY1_selection.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Fig 3: COMP1 ──────────────────────────────────────────────────────────────

def fig3_comp1():
    bm = load_bm("COMP1")
    th = load_thymus("COMP1")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1]})
    fig.suptitle("Fig 3 — COMP1: BM → Thymus (direct coupling, zero lag)", fontweight="bold")

    # BM export flux
    ax1.plot(bm["days"], bm["flux"], color=COLORS["flux"], lw=1.8, label="BM export (DN1)")
    _shade(ax1, bm["days"], bm["flux_lo"], bm["flux_hi"], COLORS["flux"])
    ax1.set_ylabel("cells · day⁻¹")
    ax1.set_title("Progenitor export from bone marrow", fontsize=9)
    ax1.legend(loc="upper right")
    ax1.grid(True, ls=":", alpha=0.4)

    # Thymus export
    ax2.plot(th["days"], th["export"], color=COLORS["CD4SP"], lw=1.8, label="Thymus naïve T export")
    _shade(ax2, th["days"], th["export_lo"], th["export_hi"], COLORS["CD4SP"])
    ax2.set_xlabel("Simulation time (days)")
    ax2.set_ylabel("cells · day⁻¹")
    ax2.set_title("Naïve T cell export from thymus", fontsize=9)
    ax2.legend(loc="upper right")
    ax2.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    out = FIGURES / "fig3_COMP1_direct.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Fig 4: COMP2 ──────────────────────────────────────────────────────────────

def fig4_comp2():
    bm = load_bm("COMP2")
    th = load_thymus("COMP2")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1]})
    fig.suptitle("Fig 4 — COMP2: BM → Blood Transit → Thymus (biologically realistic lag)",
                 fontweight="bold")

    # BM export
    ax1.plot(bm["days"], bm["flux"], color=COLORS["flux"], lw=1.8, label="BM export (DN1)")
    _shade(ax1, bm["days"], bm["flux_lo"], bm["flux_hi"], COLORS["flux"])
    ax1.set_ylabel("cells · day⁻¹")
    ax1.set_title("Progenitor export from bone marrow", fontsize=9)
    ax1.legend(loc="upper right")
    ax1.grid(True, ls=":", alpha=0.4)

    # Thymus export — annotate the lag
    ax2.plot(th["days"], th["export"], color=COLORS["CD4SP"], lw=1.8, label="Thymus naïve T export")
    _shade(ax2, th["days"], th["export_lo"], th["export_hi"], COLORS["CD4SP"])

    # Annotate blood transit lag (~7 days)
    lag_days = 601933.2 / 86400
    ax2.axvline(lag_days, color="gray", ls="--", lw=1, alpha=0.7)
    ax2.text(lag_days + 0.3, ax2.get_ylim()[1] * 0.05 if ax2.get_ylim()[1] > 0 else 0.01,
             f"Transit lag\n≈{lag_days:.1f} d", fontsize=7, color="gray", va="bottom")

    ax2.set_xlabel("Simulation time (days)")
    ax2.set_ylabel("cells · day⁻¹")
    ax2.set_title("Naïve T cell export from thymus (delayed by blood transit)", fontsize=9)
    ax2.legend(loc="upper right")
    ax2.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    out = FIGURES / "fig4_COMP2_transfer.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Fig 5: COMP3 ──────────────────────────────────────────────────────────────

def fig5_comp3():
    bm = load_bm("COMP3")
    th = load_thymus("COMP3")
    pn = load_periph("COMP3")

    fig, axes = plt.subplots(3, 1, figsize=(7, 7.5), sharex=True,
                              gridspec_kw={"height_ratios": [1, 1, 1.5]})
    fig.suptitle("Fig 5 — COMP3: Full Immune Ontogeny Pipeline\n"
                 "BM → Blood Transit → Thymus → Peripheral LN", fontweight="bold")

    ax1, ax2, ax3 = axes

    # BM export
    ax1.plot(bm["days"], bm["flux"], color=COLORS["flux"], lw=1.8, label="BM progenitor export")
    _shade(ax1, bm["days"], bm["flux_lo"], bm["flux_hi"], COLORS["flux"])
    ax1.set_ylabel("cells · day⁻¹")
    ax1.set_title("Bone marrow → thymus progenitor flux", fontsize=9)
    ax1.legend(loc="upper right")
    ax1.grid(True, ls=":", alpha=0.4)

    # Transit lag annotation
    lag_days = 601933.2 / 86400
    ax1.axvline(lag_days, color="gray", ls="--", lw=0.8, alpha=0.6)

    # Thymus export
    ax2.plot(th["days"], th["export"], color=COLORS["CD4SP"], lw=1.8, label="Thymus export")
    _shade(ax2, th["days"], th["export_lo"], th["export_hi"], COLORS["CD4SP"])
    ax2.set_ylabel("cells · day⁻¹")
    ax2.set_title("Thymus → periphery naïve T export", fontsize=9)
    ax2.legend(loc="upper right")
    ax2.grid(True, ls=":", alpha=0.4)

    # Peripheral LN pools — scale to thousands
    cd4_k = pn["cd4"] / 1000
    cd8_k = pn["cd8"] / 1000
    ax3.plot(pn["days"], cd4_k, color=COLORS["CD4"], lw=1.8, label="Naïve CD4+ T")
    ax3.plot(pn["days"], cd8_k, color=COLORS["CD8"], lw=1.8, label="Naïve CD8+ T")
    ax3.fill_between(pn["days"], cd4_k, cd8_k, alpha=0.1, color="purple")
    ax3.set_xlabel("Simulation time (days)")
    ax3.set_ylabel("Cell count (×10³)")
    ax3.set_title("Peripheral naïve T cell pools (Borghans-De Boer homeostasis)", fontsize=9)
    ax3.legend(loc="upper right")
    ax3.grid(True, ls=":", alpha=0.4)

    # Annotate homing lag (2 days after first thymus export)
    homing_lag_days = 2.0
    first_export_day = th["days"][th["export"] > 0][0] if np.any(th["export"] > 0) else 0
    if first_export_day > 0:
        ax2.axvline(first_export_day, color="purple", ls=":", lw=0.8, alpha=0.5)
        ax3.axvline(first_export_day + homing_lag_days, color="purple", ls=":", lw=0.8, alpha=0.5)

    fig.tight_layout()
    out = FIGURES / "fig5_COMP3_full_graph.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Generating figures…")
    fig1_bm1()
    fig2_thy1()
    fig3_comp1()
    fig4_comp2()
    fig5_comp3()
    print("Done.")
