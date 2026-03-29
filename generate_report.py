#!/usr/bin/env python3
"""
OISA Simulation Report Generator
==================================
Reads ISSL checkpoints from results/, (optionally) regenerates figures via
plot_results.py, then assembles a self-contained PDF report with:
  - Cover page + framework overview
  - Per-simulation section: architecture diagram, description, metrics table,
    simulation figure
  - Comparative summary

Usage
-----
    python generate_report.py                          # figs + PDF
    python generate_report.py --no-figs               # skip figure regen
    python generate_report.py --output report.pdf
"""

from __future__ import annotations

import argparse
import io
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

# ── ReportLab ─────────────────────────────────────────────────────────────────
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    HRFlowable,
    Image,
    KeepTogether,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

# ── matplotlib ────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO      = Path(__file__).parent
RESULTS   = REPO / "results"
FIGS_DIR  = RESULTS / "figures"
OUTPUT    = RESULTS / "OISA_simulation_report.pdf"

# ── Colour palette ────────────────────────────────────────────────────────────
C_NAVY    = colors.HexColor("#1A3A5C")
C_BLUE    = colors.HexColor("#2E75B6")
C_ORANGE  = colors.HexColor("#C05621")
C_GREEN   = colors.HexColor("#276749")
C_PURPLE  = colors.HexColor("#553C9A")
C_LGRAY   = colors.HexColor("#F5F7FA")
C_TBLHDR  = colors.HexColor("#2E75B6")
C_TBLALT  = colors.HexColor("#EBF3FB")
C_TEXT    = colors.HexColor("#2D3748")
C_MUTED   = colors.HexColor("#718096")
WHITE     = colors.white
BLACK     = colors.black

# Diagram colours per formalism
DGM_COLORS = {
    "ODE":      "#2E75B6",
    "ABM":      "#C05621",
    "transfer": "#276749",
}

# ── Entity → human label mapping ──────────────────────────────────────────────
ENTITY_LABELS: dict[str, str] = {
    "CL:0000037": "HSC",
    "CL:0000837": "MPP",
    "CL:0000838": "LMPP",
    "CL:0000051": "CLP",
    "CL:0002420": "DN1 / ETP",
    "CL:0000893": "DP thymocyte",
    "CL:0000624": "CD4⁺ T",
    "CL:0000625": "CD8⁺ T",
    "CL:0000898": "Naïve T (total)",
}

# ── Per-scenario metadata ─────────────────────────────────────────────────────
SCENARIOS: list[dict[str, Any]] = [
    {
        "key":      "BM1_baseline",
        "title":    "Simulation 1 — Bone Marrow Haematopoiesis Baseline",
        "subtitle": "§5.1 · Standalone ODE · 30 days · Δt = 6 h",
        "figure":   "fig1_BM1_baseline.png",
        "models":   ["bm_haematopoiesis"],
        "nodes":    [("Bone Marrow\nODE", "Δt = 6 h", "ODE")],
        "edges":    [],
        "description": (
            "The bone marrow model simulates progenitor production through a "
            "five-compartment ODE cascade: <b>HSC → MPP → LMPP → CLP → DN1</b>. "
            "HSC renewal follows logistic dynamics (K<sub>niche</sub> = 11 000). "
            "MPP operates in a <i>near-balanced</i> proliferation regime "
            "(λ<sub>MPP</sub> = d<sub>MPP_LMPP</sub> + d<sub>MPP_death</sub> − r<sub>MPP</sub> = 0.001 day⁻¹) "
            "that amplifies MPP count relative to input while maintaining stability. "
            "Initial conditions are computed analytically from the parameter set "
            "(v7) so the model starts exactly at steady state, producing a flat "
            "export flux from t = 0. "
            "At equilibrium the DN1 export flux is <b>~59.6 cells·day⁻¹</b>, "
            "consistent with the murine literature range of 10–100 cells·day⁻¹ "
            "(Bhandoola et al. 2007). "
            "CI-95 bands are computed by Monte Carlo parameter sampling "
            "(2 000 draws, ~51 % stable), conditioning on λ<sub>MPP</sub> &gt; 0; "
            "the resulting flux CI [1.4, 90.9] spans the full literature range, "
            "reflecting genuine parametric uncertainty of the near-balanced regime."
        ),
    },
    {
        "key":      "THY1_baseline",
        "title":    "Simulation 2 — Thymus Selection Baseline",
        "subtitle": "§5.2 · Standalone ABM · 30 days · Δt = 24 h",
        "figure":   "fig2_THY1_baseline.png",
        "models":   ["thymus_selection"],
        "nodes":    [("Thymus\nABM", "Δt = 24 h", "ABM")],
        "edges":    [],
        "description": (
            "The thymus model is an <b>agent-based simulation</b> of thymocyte "
            "development and TCR selection. Each agent represents "
            "<b>300 000 real thymocytes</b> (scaling factor, Scollay &amp; Godfrey 1995). "
            "Agents progress through: DN1 import → DN→DP transition (~20 substeps "
            "at 1 h each) → cortical positive selection (TCR affinity threshold) → "
            "medullary dwell (<i>medullary_dwell_steps</i> = 72 substeps = 3 days) "
            "→ export as naïve T cells. "
            "In standalone mode, progenitor import is set to "
            "<b>baseline_import = 60 cells·day⁻¹</b> "
            "(calibrated to the BM v7 equilibrium flux). "
            "12 independent realisations are run per checkpoint; CI-95 is computed "
            "empirically across realisations. "
            "After the initial maturation lag (~4–5 days), the model exports "
            "~1 × 10⁶ cells·day⁻¹ (scaled), consistent with murine thymic output. "
            "Naïve T export is available in <code>export_signals[0].biological_flux_per_day</code>."
        ),
    },
    {
        "key":      "COMP1_direct",
        "title":    "Simulation 3 — BM → Thymus Direct Coupling",
        "subtitle": "§5.3a · ODE → ABM · 30 days",
        "figure":   "fig3_COMP1_direct.png",
        "models":   ["bm_haematopoiesis", "thymus_selection"],
        "nodes":    [
            ("Bone Marrow\nODE", "Δt = 6 h",  "ODE"),
            ("Thymus\nABM",      "Δt = 24 h", "ABM"),
        ],
        "edges":    [(0, 1, "progenitor_export\n(direct, lag = 0)")],
        "description": (
            "The first coupled scenario connects the BM export signal directly "
            "to the thymus import, with <b>no transit delay</b>. "
            "The thymus <code>baseline_import</code> override is set to 0 so "
            "that all progenitor input comes from the live BM signal. "
            "The BM step runs every 6 h; for each 24 h thymus checkpoint the "
            "most recent BM signal is used as import. "
            "This scenario establishes the baseline coupling behaviour and shows "
            "how the BM-calibrated flux (~59.6 cells·day⁻¹) feeds thymic "
            "population dynamics without physiological transport delay. "
            "The slightly lower import compared to the standalone baseline (60 cells·day⁻¹) "
            "produces marginally different thymocyte pool kinetics, "
            "visible in the DN1 and DP panels."
        ),
    },
    {
        "key":      "COMP2_transfer",
        "title":    "Simulation 4 — BM → Blood Transit → Thymus",
        "subtitle": "§5.3b · ODE → transfer → ABM · 30 days",
        "figure":   "fig4_COMP2_transfer.png",
        "models":   ["bm_haematopoiesis", "blood_transit", "thymus_selection"],
        "nodes":    [
            ("Bone Marrow\nODE",      "Δt = 6 h",       "ODE"),
            ("Blood Transit\ntransfer", "on-demand",     "transfer"),
            ("Thymus\nABM",           "Δt = 24 h",      "ABM"),
        ],
        "edges":    [
            (0, 1, "progenitor_export"),
            (1, 2, "delivery + lag ~4 d\n(stop_fraction = 0.82)"),
        ],
        "description": (
            "The blood transit model inserts a <b>physiologically realistic "
            "~4-day delay</b> between BM export and thymic arrival, "
            "modelling the time progenitors spend circulating in the bloodstream "
            "before homing to the thymus. "
            "The transfer ODE is stateless and invoked on-demand: it receives "
            "the BM signal, computes a retention fraction "
            "(<i>stop_fraction</i> = 0.82) and a lag time in seconds, "
            "and enqueues the attenuated signal for thymus delivery. "
            "The orchestrator (or in-process runner) holds each signal in a "
            "pending queue and delivers it once <code>sim_time ≥ deliver_at</code>. "
            "Comparing COMP2 to COMP1 isolates the effect of the transit delay "
            "on thymic engraftment kinetics: thymocyte pools reach steady state "
            "~4 days later, and the initial DN1 pool grows more gradually."
        ),
    },
    {
        "key":      "COMP3_full_graph",
        "title":    "Simulation 5 — Full Immune Ontogeny Pipeline",
        "subtitle": "§5.4 · ODE → transfer → ABM → ODE · 30 days",
        "figure":   "fig5_COMP3_full_graph.png",
        "models":   ["bm_haematopoiesis", "blood_transit", "thymus_selection", "peripheral_ln"],
        "nodes":    [
            ("Bone Marrow\nODE",      "Δt = 6 h",    "ODE"),
            ("Blood Transit\ntransfer", "on-demand",  "transfer"),
            ("Thymus\nABM",           "Δt = 24 h",   "ABM"),
            ("Peripheral LN\nODE",    "Δt = 12 h",   "ODE"),
        ],
        "edges":    [
            (0, 1, "progenitor_export"),
            (1, 2, "+lag ~4 d"),
            (2, 3, "naive_T_export\n+lag 2 d (homing)"),
        ],
        "description": (
            "The complete immune ontogeny pipeline chains all four models. "
            "Naive T cells exported by the thymus are routed to the "
            "<b>peripheral lymph node (PLN)</b> model with a "
            "<b>2-day homing lag</b> (constant, declared in the config graph). "
            "The PLN implements Borghans–De Boer homeostatic dynamics for "
            "CD4⁺ and CD8⁺ naïve T pools: "
            "d[CD4]/dt = import<sub>CD4</sub> + ρ₄(S₄ − CD4) − d<sub>total</sub>·CD4. "
            "Import is routed via <code>lymph_node_fraction = 0.0013</code> "
            "(fraction of total thymic output assigned to this compartment). "
            "Without thymic input the Borghans equilibrium converges to S/2 "
            "(τ ≈ 167 days); with correct import the pools are maintained near "
            "their set points (S₄ = 200 000, S₈ = 100 000). "
            "The CD4/CD8 ratio evolves from 2.00 toward ~1.92 "
            "as thymic export builds up, reflecting the cd4_fraction = 0.65 "
            "split from the thymus."
        ),
    },
]

# ═══════════════════════════════════════════════════════════════════════════════
# ISSL data loading
# ═══════════════════════════════════════════════════════════════════════════════

def _last_checkpoint(issl_model_dir: Path) -> dict | None:
    files = sorted(issl_model_dir.glob("checkpoint_*.json"),
                   key=lambda f: int(f.stem.split("_")[1]))
    return json.loads(files[-1].read_text()) if files else None


def load_scenario(key: str) -> dict[str, dict]:
    """Return {model_id: last_checkpoint_dict} for a scenario."""
    base = RESULTS / key / key / "issl"
    result = {}
    if not base.exists():
        return result
    for model_dir in base.iterdir():
        ck = _last_checkpoint(model_dir)
        if ck:
            result[model_dir.name] = ck
    return result


def fmt_num(v: float, digits: int = 1) -> str:
    if abs(v) >= 1_000_000:
        return f"{v/1_000_000:.{digits}f} M"
    if abs(v) >= 1_000:
        return f"{v:,.{digits}f}"
    return f"{v:.{digits}f}"


def fmt_ci(ci: list[float]) -> str:
    if not ci or len(ci) < 2:
        return "—"
    lo, hi = ci[0], ci[1]
    if lo < 0 or abs(lo) > 1e8:
        return "N/A *"
    return f"[{fmt_num(lo, 0)}, {fmt_num(hi, 0)}]"


def extract_metrics(scenario_data: dict[str, dict]) -> dict[str, list[list[str]]]:
    """Return per-model metric rows: [[name, value, ci, unit], ...]."""
    rows: dict[str, list[list[str]]] = {}

    # ── BM ─────────────────────────────────────────────────────────────────────
    bm = scenario_data.get("bm_haematopoiesis")
    if bm:
        cs   = {e["entity_id"]: e for e in bm.get("continuous_state", [])}
        sig  = bm.get("export_signals", [{}])[0]
        t_d  = bm["envelope"]["sim_time_s"] / 86400
        r    = [["Compartment", f"Day {t_d:.0f} (cells)", "CI-95 (MC)", ""]]
        for eid, label in [
            ("CL:0000037", "HSC"),
            ("CL:0000837", "MPP"),
            ("CL:0000838", "LMPP"),
            ("CL:0000051", "CLP"),
            ("CL:0002420", "DN1 / ETP"),
        ]:
            e = cs.get(eid, {})
            r.append([label, fmt_num(e.get("count", 0)),
                      fmt_ci(e.get("ci_95", [])), "cells"])
        r.append(["Export flux",
                  f"{sig.get('flux', 0):.2f}",
                  fmt_ci(sig.get("ci_95", [])),
                  "cells·day⁻¹"])
        rows["Bone Marrow"] = r

    # ── Thymus ─────────────────────────────────────────────────────────────────
    thy = scenario_data.get("thymus_selection")
    if thy:
        cs   = {e["entity_id"]: e for e in thy.get("continuous_state", [])}
        sig  = thy.get("export_signals", [{}])[0]
        t_d  = thy["envelope"]["sim_time_s"] / 86400
        sf   = sig.get("scaling_factor", 1)
        r    = [["Population", f"Day {t_d:.0f} (cells, scaled)", "CI-95", ""]]
        for eid, label in [
            ("CL:0002420", "DN1 / ETP"),
            ("CL:0000893", "DP thymocyte"),
            ("CL:0000624", "CD4⁺ SP"),
            ("CL:0000625", "CD8⁺ SP"),
        ]:
            e = cs.get(eid, {})
            r.append([label, fmt_num(e.get("count", 0)),
                      fmt_ci(e.get("ci_95", [])), "cells"])
        bio = sig.get("biological_flux_per_day", sig.get("flux", 0) * sf)
        r.append(["Naïve T export",
                  fmt_num(bio),
                  fmt_ci(sig.get("ci_95", [])),
                  "cells·day⁻¹"])
        rows["Thymus"] = r

    # ── Blood Transit ──────────────────────────────────────────────────────────
    bt = scenario_data.get("blood_transit")
    if bt:
        cs   = {e["entity_id"]: e for e in bt.get("continuous_state", [])}
        sig  = bt.get("export_signals", [{}])[0]
        lag  = sig.get("lag_s", 0) or 0
        r = [["Metric", "Value", "CI-95", ""]]
        circ = cs.get("CL:0002420", {})
        r.append(["Circulating progenitors",
                  fmt_num(circ.get("count", 0)),
                  fmt_ci(circ.get("ci_95", [])), "cells"])
        r.append(["Delivered flux",
                  f"{sig.get('flux', 0):.2f}",
                  fmt_ci(sig.get("ci_95", [])), "cells·day⁻¹"])
        r.append(["Transit lag", f"{lag/3600:.1f}", "—", "hours"])
        rows["Blood Transit"] = r

    # ── PLN ────────────────────────────────────────────────────────────────────
    pln = scenario_data.get("peripheral_ln")
    if pln:
        cs  = {e["entity_id"]: e for e in pln.get("continuous_state", [])}
        t_d = pln["envelope"]["sim_time_s"] / 86400
        cd4 = cs.get("CL:0000624", {}).get("count", 0)
        cd8 = cs.get("CL:0000625", {}).get("count", 0)
        r = [["Pool", f"Day {t_d:.0f} (cells)", "Set point", ""]]
        r.append(["CD4⁺ naïve T", fmt_num(cd4, 0), "200 000", "cells"])
        r.append(["CD8⁺ naïve T", fmt_num(cd8, 0), "100 000", "cells"])
        ratio = cd4 / cd8 if cd8 > 0 else 0
        r.append(["CD4/CD8 ratio", f"{ratio:.3f}", "2.000 (initial)", "—"])
        rows["Peripheral LN"] = r

    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline diagram generator
# ═══════════════════════════════════════════════════════════════════════════════

def draw_pipeline(nodes: list[tuple[str, str, str]],
                  edges: list[tuple[int, int, str]],
                  width_in: float = 13.0) -> bytes:
    """
    Draw a pipeline diagram as PNG bytes.
    nodes: [(main_label, sub_label, formalism), ...]
    edges: [(from_idx, to_idx, edge_label), ...]
    """
    n     = len(nodes)
    h_in  = 2.0 if n <= 3 else 2.2
    fig, ax = plt.subplots(figsize=(width_in, h_in), facecolor="white")
    ax.set_xlim(-0.3, n - 1 + 0.3 + (n - 1) * 0.4)
    ax.set_ylim(-0.6, 0.7)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.08)

    BOX_W, BOX_H = 1.5, 0.5
    SPACING      = 2.2 if n <= 3 else 1.9

    positions = [i * SPACING for i in range(n)]

    for i, (mlabel, slabel, formalism) in enumerate(nodes):
        xc = positions[i]
        clr = DGM_COLORS.get(formalism, "#888888")
        box = FancyBboxPatch(
            (xc - BOX_W / 2, -BOX_H / 2), BOX_W, BOX_H,
            boxstyle="round,pad=0.06",
            facecolor=clr, edgecolor="white", linewidth=1.5, alpha=0.93,
            zorder=3,
        )
        ax.add_patch(box)
        lines = mlabel.split("\n")
        ax.text(xc, 0.10 if len(lines) > 1 else 0.02, lines[0],
                ha="center", va="center",
                fontsize=10, fontweight="bold", color="white", zorder=4)
        if len(lines) > 1:
            ax.text(xc, -0.13, lines[1],
                    ha="center", va="center",
                    fontsize=8, color="white", alpha=0.9, zorder=4)
        ax.text(xc, -BOX_H / 2 - 0.12, slabel,
                ha="center", va="top",
                fontsize=7.5, color="#4A5568", style="italic", zorder=4)

    for (fi, ti, elabel) in edges:
        x0 = positions[fi] + BOX_W / 2
        x1 = positions[ti] - BOX_W / 2
        mx = (x0 + x1) / 2
        ax.annotate(
            "", xy=(x1, 0), xytext=(x0, 0),
            arrowprops=dict(
                arrowstyle="->", color="#4A5568",
                lw=1.8, mutation_scale=16,
            ),
            zorder=2,
        )
        for j, line in enumerate(elabel.split("\n")):
            ax.text(mx, 0.30 - j * 0.17, line,
                    ha="center", va="bottom",
                    fontsize=7.5, color="#4A5568", zorder=4)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def draw_full_overview() -> bytes:
    """Full 4-model pipeline diagram for the overview page."""
    nodes = [
        ("Bone Marrow",    "ODE · Δt = 6 h",     "ODE"),
        ("Blood Transit",  "transfer · on-demand","transfer"),
        ("Thymus",         "ABM · Δt = 24 h",     "ABM"),
        ("Peripheral LN",  "ODE · Δt = 12 h",     "ODE"),
    ]
    edges = [
        (0, 1, "progenitor_export"),
        (1, 2, "+lag ~4 d"),
        (2, 3, "naive_T_export\n+lag 2 d"),
    ]
    return draw_pipeline(nodes, edges, width_in=13.0)


# ═══════════════════════════════════════════════════════════════════════════════
# ReportLab styles
# ═══════════════════════════════════════════════════════════════════════════════

def make_styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    S: dict[str, ParagraphStyle] = {}

    def s(name, **kw) -> ParagraphStyle:
        p = ParagraphStyle(name, parent=base["Normal"], **kw)
        S[name] = p
        return p

    s("title_cover",    fontSize=28, leading=34, textColor=WHITE,
      fontName="Helvetica-Bold", alignment=TA_LEFT, spaceAfter=6)
    s("subtitle_cover", fontSize=13, leading=16, textColor=colors.HexColor("#BEE3F8"),
      fontName="Helvetica", alignment=TA_LEFT, spaceAfter=4)
    s("meta_cover",     fontSize=10, leading=13, textColor=colors.HexColor("#90CDF4"),
      fontName="Helvetica", alignment=TA_LEFT)

    s("h1",  fontSize=17, leading=22, textColor=C_NAVY,
      fontName="Helvetica-Bold", spaceBefore=14, spaceAfter=6)
    s("h1sub", fontSize=10, leading=13, textColor=C_MUTED,
      fontName="Helvetica", spaceAfter=10)
    s("h2",  fontSize=13, leading=17, textColor=C_BLUE,
      fontName="Helvetica-Bold", spaceBefore=10, spaceAfter=5)
    s("body", fontSize=9.5, leading=14, textColor=C_TEXT,
      fontName="Helvetica", spaceAfter=6, alignment=TA_LEFT)
    s("caption", fontSize=8.5, leading=11, textColor=C_MUTED,
      fontName="Helvetica-Oblique", spaceAfter=4, alignment=TA_CENTER)
    s("hdr",  fontSize=8, leading=10, textColor=C_MUTED,
      fontName="Helvetica", alignment=TA_RIGHT)
    s("ftr",  fontSize=8, leading=10, textColor=C_MUTED,
      fontName="Helvetica", alignment=TA_CENTER)
    s("tbl_hdr",  fontSize=8.5, leading=11, textColor=WHITE,
      fontName="Helvetica-Bold", alignment=TA_LEFT)
    s("tbl_cell", fontSize=8.5, leading=11, textColor=C_TEXT,
      fontName="Helvetica", alignment=TA_LEFT)
    s("note", fontSize=8, leading=11, textColor=C_MUTED,
      fontName="Helvetica-Oblique", spaceAfter=4)

    return S


def _tbl_style(n_rows: int) -> TableStyle:
    cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), C_TBLHDR),
        ("TEXTCOLOR",  (0, 0), (-1, 0), WHITE),
        ("FONTNAME",   (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",   (0, 0), (-1, -1), 8.5),
        ("ROWBACKGROUND", (0, 1), (-1, -1),
         [C_LGRAY if i % 2 == 0 else WHITE for i in range(n_rows - 1)]),
        ("GRID",       (0, 0), (-1, -1), 0.25, colors.HexColor("#CBD5E0")),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING",   (0, 0), (-1, -1), 7),
        ("RIGHTPADDING",  (0, 0), (-1, -1), 7),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
    ]
    return TableStyle(cmds)


# ═══════════════════════════════════════════════════════════════════════════════
# Canvas callbacks (header / footer)
# ═══════════════════════════════════════════════════════════════════════════════

PAGE_W, PAGE_H = A4

def _on_cover(canvas, doc):
    canvas.saveState()
    # Deep navy gradient background
    canvas.setFillColor(C_NAVY)
    canvas.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    # Accent stripe
    canvas.setFillColor(C_BLUE)
    canvas.rect(0, PAGE_H * 0.38, PAGE_W, PAGE_H * 0.62, fill=1, stroke=0)
    # Bottom stripe
    canvas.setFillColor(colors.HexColor("#0D2137"))
    canvas.rect(0, 0, PAGE_W, PAGE_H * 0.12, fill=1, stroke=0)
    canvas.restoreState()


def _on_content(canvas, doc):
    canvas.saveState()
    # Header line
    canvas.setStrokeColor(colors.HexColor("#E2E8F0"))
    canvas.setLineWidth(0.5)
    canvas.line(2 * cm, PAGE_H - 1.5 * cm, PAGE_W - 2 * cm, PAGE_H - 1.5 * cm)
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(C_MUTED)
    canvas.drawRightString(PAGE_W - 2 * cm, PAGE_H - 1.3 * cm,
                           "OISA — Immune Ontogeny Simulation Report")
    # Footer
    canvas.line(2 * cm, 1.5 * cm, PAGE_W - 2 * cm, 1.5 * cm)
    canvas.drawCentredString(PAGE_W / 2, 1.15 * cm, str(doc.page - 1))
    canvas.restoreState()


# ═══════════════════════════════════════════════════════════════════════════════
# Flowable helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _img_from_bytes(data: bytes, width_pt: float) -> Image:
    buf = io.BytesIO(data)
    img = Image(buf)
    aspect = img.imageHeight / img.imageWidth
    img._restrictSize(width_pt, width_pt * aspect)
    return img


def _fig_image(filename: str, width_pt: float) -> Image | None:
    p = FIGS_DIR / filename
    if not p.exists():
        return None
    img = Image(str(p))
    aspect = img.imageHeight / img.imageWidth
    img._restrictSize(width_pt, width_pt * aspect)
    return img


def _metrics_table(rows: list[list[str]], S: dict) -> Table:
    styled = []
    for i, row in enumerate(rows):
        if i == 0:
            styled.append([Paragraph(c, S["tbl_hdr"]) for c in row])
        else:
            styled.append([Paragraph(c, S["tbl_cell"]) for c in row])
    col_w = [5.5 * cm, 4.0 * cm, 4.5 * cm, 2.5 * cm]
    t = Table(styled, colWidths=col_w)
    t.setStyle(_tbl_style(len(rows)))
    return t


def _hr(S) -> HRFlowable:
    return HRFlowable(width="100%", thickness=0.5,
                      color=colors.HexColor("#CBD5E0"),
                      spaceAfter=6, spaceBefore=6)


def _legend_table(S: dict) -> Table:
    items = [
        ("■", DGM_COLORS["ODE"],      "ODE — Ordinary Differential Equation"),
        ("■", DGM_COLORS["ABM"],      "ABM — Agent-Based Model"),
        ("■", DGM_COLORS["transfer"], "Transfer — stateless delay model"),
    ]
    row = []
    for sym, clr, label in items:
        row.append(Paragraph(
            f'<font color="{clr}"><b>{sym}</b></font>  {label}',
            S["note"],
        ))
    t = Table([row], colWidths=[5.5 * cm, 5.5 * cm, 5.5 * cm])
    t.setStyle(TableStyle([
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    return t


# ═══════════════════════════════════════════════════════════════════════════════
# Cover page content
# ═══════════════════════════════════════════════════════════════════════════════

def cover_flowables(S: dict) -> list:
    TW = PAGE_W - 4 * cm   # text width
    story = []
    # Vertical push to lower third
    story.append(Spacer(TW, 9.5 * cm))
    story.append(Paragraph("OISA", S["title_cover"]))
    story.append(Paragraph("Immune Ontogeny Simulation", S["title_cover"]))
    story.append(Spacer(TW, 0.4 * cm))
    story.append(Paragraph("Simulation Report", S["subtitle_cover"]))
    story.append(Spacer(TW, 0.3 * cm))
    story.append(Paragraph(
        "Five coupled scenarios · BM → Blood Transit → Thymus → Peripheral LN · 30 days",
        S["subtitle_cover"],
    ))
    story.append(Spacer(TW, 1.2 * cm))
    from datetime import date
    story.append(Paragraph(f"Generated: {date.today().isoformat()}", S["meta_cover"]))
    story.append(Paragraph(
        "Models: BM v7 · Thymus v3 · Blood Transit · PLN", S["meta_cover"]
    ))
    return story


# ═══════════════════════════════════════════════════════════════════════════════
# Overview page
# ═══════════════════════════════════════════════════════════════════════════════

def overview_flowables(S: dict) -> list:
    TW = PAGE_W - 4 * cm
    story: list = []
    story.append(Paragraph("OISA Framework Overview", S["h1"]))
    story.append(Paragraph(
        "Orchestrated Immune Simulation Architecture", S["h1sub"],
    ))
    story.append(Paragraph(
        "OISA couples heterogeneous computational models — ODEs and an ABM — "
        "into a directed graph where each node is a <b>domain model</b> and "
        "each edge is an <b>ISSL signal</b> (Internal Simulation State Log). "
        "The orchestrator advances models on their own time steps, routes "
        "export signals with optional lags, and reconciles the stochastic ABM "
        "output with deterministic ODE inputs by using the distributional mean. "
        "This report covers five progressively coupled scenarios, from single-model "
        "baselines to the full four-model immune ontogeny pipeline.",
        S["body"],
    ))
    story.append(Spacer(TW, 0.3 * cm))
    story.append(Paragraph("Full pipeline — COMP3", S["h2"]))
    dgm = draw_full_overview()
    story.append(_img_from_bytes(dgm, TW))
    story.append(Spacer(TW, 0.2 * cm))
    story.append(_legend_table(S))
    story.append(Spacer(TW, 0.5 * cm))
    story.append(_hr(S))
    story.append(Paragraph("Simulation scenarios", S["h2"]))

    tbl_data = [
        [Paragraph(h, S["tbl_hdr"]) for h in
         ["#", "Scenario", "Models", "Duration", "Δt (fastest)"]],
        *[
            [Paragraph(c, S["tbl_cell"]) for c in row]
            for row in [
                ["1", "BM1 Baseline",      "BM",                              "30 d", "6 h"],
                ["2", "THY1 Baseline",     "Thymus",                          "30 d", "24 h"],
                ["3", "COMP1 Direct",      "BM → Thymus",                     "30 d", "6 h"],
                ["4", "COMP2 Transfer",    "BM → Blood Transit → Thymus",     "30 d", "6 h"],
                ["5", "COMP3 Full Graph",  "BM → Transit → Thymus → PLN",     "30 d", "6 h"],
            ]
        ],
    ]
    col_w = [0.8*cm, 4.2*cm, 7.0*cm, 2.0*cm, 2.5*cm]
    tbl = Table(tbl_data, colWidths=col_w)
    tbl.setStyle(_tbl_style(6))
    story.append(tbl)
    return story


# ═══════════════════════════════════════════════════════════════════════════════
# Per-scenario section builder
# ═══════════════════════════════════════════════════════════════════════════════

def scenario_flowables(meta: dict, S: dict) -> list:
    TW = PAGE_W - 4 * cm
    story: list = []

    story.append(Paragraph(meta["title"],    S["h1"]))
    story.append(Paragraph(meta["subtitle"], S["h1sub"]))

    # Architecture diagram
    story.append(Paragraph("Pipeline topology", S["h2"]))
    dgm_bytes = draw_pipeline(meta["nodes"], meta["edges"])
    diag_w    = min(TW, len(meta["nodes"]) * 4.0 * cm + 1.0 * cm)
    story.append(KeepTogether([
        _img_from_bytes(dgm_bytes, diag_w),
        Spacer(TW, 0.1 * cm),
        _legend_table(S),
    ]))

    story.append(Spacer(TW, 0.4 * cm))

    # Description
    story.append(Paragraph("Description", S["h2"]))
    story.append(Paragraph(meta["description"], S["body"]))

    story.append(Spacer(TW, 0.4 * cm))

    # Metrics tables
    data = load_scenario(meta["key"])
    if data:
        metrics = extract_metrics(data)
        if metrics:
            story.append(Paragraph("Key metrics — day 30", S["h2"]))
            for model_name, rows in metrics.items():
                story.append(Paragraph(model_name, S["note"]))
                story.append(_metrics_table(rows, S))
                story.append(Spacer(TW, 0.2 * cm))
            if any("N/A" in str(r) for rows in metrics.values() for r in rows):
                story.append(Paragraph(
                    "* CI-95 'N/A' indicates the analytical formula is unreliable "
                    "for this compartment. Monte Carlo CI is not yet implemented "
                    "for this model.",
                    S["note"],
                ))

    story.append(Spacer(TW, 0.4 * cm))

    # Simulation figure
    fig = _fig_image(meta["figure"], TW)
    if fig:
        story.append(Paragraph("Simulation output", S["h2"]))
        story.append(fig)
        story.append(Paragraph(
            f"Figure — {meta['title']}. "
            "Compartment trajectories over 30 days; shaded bands = CI-95.",
            S["caption"],
        ))

    return story


# ═══════════════════════════════════════════════════════════════════════════════
# Comparative section
# ═══════════════════════════════════════════════════════════════════════════════

def comparative_flowables(S: dict) -> list:
    TW = PAGE_W - 4 * cm
    story: list = []
    story.append(Paragraph("Comparative Summary", S["h1"]))
    story.append(Paragraph(
        "Cross-scenario overview — DN1 export flux and naïve T output at day 30",
        S["h1sub"],
    ))
    story.append(Paragraph(
        "The table below compares the key output metrics across the five scenarios. "
        "The progressive addition of transit delay (COMP2) and peripheral homeostasis "
        "(COMP3) reveals the cumulative impact of each biological layer on naive T-cell "
        "pool dynamics.",
        S["body"],
    ))

    # Build comparison table from last ISSL checkpoints
    rows = [
        [Paragraph(h, S["tbl_hdr"]) for h in
         ["Scenario", "BM flux\n(cells·day⁻¹)", "Thymus export\n(M cells·day⁻¹)",
          "Transit lag\n(h)", "CD4 pool", "CD4/CD8"]],
    ]
    for meta in SCENARIOS:
        data = load_scenario(meta["key"])
        bm_f   = "—"; thy_f = "—"; lag = "—"; cd4 = "—"; ratio = "—"
        if "bm_haematopoiesis" in data:
            sig  = data["bm_haematopoiesis"].get("export_signals", [{}])[0]
            bm_f = f"{sig.get('flux', 0):.1f}"
        if "thymus_selection" in data:
            sig    = data["thymus_selection"].get("export_signals", [{}])[0]
            bio    = sig.get("biological_flux_per_day") or 0
            thy_f  = f"{bio/1e6:.2f}" if bio else "—"
        if "blood_transit" in data:
            sig = data["blood_transit"].get("export_signals", [{}])[0]
            l   = sig.get("lag_s", 0) or 0
            lag = f"{l/3600:.0f}"
        if "peripheral_ln" in data:
            cs  = {e["entity_id"]: e for e in
                   data["peripheral_ln"].get("continuous_state", [])}
            c4  = cs.get("CL:0000624", {}).get("count", 0)
            c8  = cs.get("CL:0000625", {}).get("count", 0)
            cd4 = fmt_num(c4, 0)
            ratio = f"{c4/c8:.3f}" if c8 > 0 else "—"
        rows.append([Paragraph(c, S["tbl_cell"]) for c in
                     [meta["key"].replace("_", " "), bm_f, thy_f, lag, cd4, ratio]])

    col_w = [3.5*cm, 3.0*cm, 3.5*cm, 2.5*cm, 3.0*cm, 2.0*cm]
    t = Table(rows, colWidths=col_w)
    t.setStyle(_tbl_style(len(rows)))
    story.append(t)
    story.append(Spacer(TW, 0.5 * cm))

    fig = _fig_image("fig6_comparative.png", TW)
    if fig:
        story.append(fig)
        story.append(Paragraph(
            "Figure 6 — Comparative trajectories. "
            "DN1 export flux (top row) and thymocyte populations (bottom row) "
            "across all five scenarios.",
            S["caption"],
        ))

    return story


# ═══════════════════════════════════════════════════════════════════════════════
# PDF assembly
# ═══════════════════════════════════════════════════════════════════════════════

def build_report(output: Path, regen_figs: bool = True) -> None:
    if regen_figs:
        print("Regenerating figures …")
        venv_py = REPO / ".venv" / "bin" / "python"
        py = str(venv_py) if venv_py.exists() else sys.executable
        res = subprocess.run(
            [py, "results/plot_results.py"], cwd=str(REPO),
            capture_output=True, text=True,
        )
        if res.returncode != 0:
            print("  plot_results.py stderr:", res.stderr[:400])
        else:
            print(res.stdout.strip())

    print(f"Building PDF → {output} …")
    S = make_styles()

    # Frames
    cover_frame   = Frame(2*cm, 2*cm, PAGE_W-4*cm, PAGE_H-4*cm, id="cover")
    content_frame = Frame(2*cm, 2.2*cm, PAGE_W-4*cm, PAGE_H-4*cm, id="content")

    doc = BaseDocTemplate(
        str(output),
        pagesize=A4,
        title="OISA Simulation Report",
        author="OISA",
    )
    doc.addPageTemplates([
        PageTemplate(id="Cover",   frames=cover_frame,   onPage=_on_cover),
        PageTemplate(id="Content", frames=content_frame, onPage=_on_content),
    ])

    story: list = []

    # ── Cover ────────────────────────────────────────────────────────────────
    story += cover_flowables(S)
    story.append(NextPageTemplate("Content"))
    story.append(PageBreak())

    # ── Overview ─────────────────────────────────────────────────────────────
    story += overview_flowables(S)
    story.append(PageBreak())

    # ── Five scenario sections ────────────────────────────────────────────────
    for i, meta in enumerate(SCENARIOS):
        story += scenario_flowables(meta, S)
        story.append(PageBreak())

    # ── Comparative ──────────────────────────────────────────────────────────
    story += comparative_flowables(S)

    doc.build(story)
    print(f"  ✓  {output}  ({output.stat().st_size // 1024} KB)")


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(description="OISA — PDF report generator")
    ap.add_argument("--no-figs",  action="store_true",
                    help="Skip figure regeneration (use existing PNGs)")
    ap.add_argument("--output",   default=str(OUTPUT),
                    help=f"Output PDF path (default: {OUTPUT})")
    args = ap.parse_args()
    build_report(Path(args.output), regen_figs=not args.no_figs)


if __name__ == "__main__":
    main()
