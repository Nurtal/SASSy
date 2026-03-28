#!/usr/bin/env python3
"""
In-process simulation runner + figure generator for all 5 OISA scenarios.

Bypasses ZMQ orchestration to run each model directly in-process.
Writes ISSL checkpoints to results/ then regenerates all figures.

Usage:
    python run_and_plot.py          # all scenarios
    python run_and_plot.py --figs   # figures only (from existing results)
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

# ── patch ZMQ + jsonschema BEFORE model imports ────────────────────────────────

_mock_zmq = MagicMock()
_mock_zmq.Context.return_value.socket.return_value = MagicMock()
sys.modules["zmq"] = _mock_zmq

import jsonschema
if not hasattr(jsonschema, "Draft202012Validator"):
    jsonschema.Draft202012Validator = jsonschema.Draft7Validator

# ── now safe to import models ─────────────────────────────────────────────────

sys.path.insert(0, str(Path(__file__).parent))

from models._base.issl_writer import write_checkpoint
from models.bm_haematopoiesis.model import BMHaematopoiesis
from models.blood_transit.model import BloodTransit
from models.peripheral_ln.model import PeripheralLN
from models.thymus_selection.model import ThymusSelection

RESULTS = Path(__file__).parent / "results"

# ── helpers ───────────────────────────────────────────────────────────────────

def _make_bm(run_issl_dir: Path) -> BMHaematopoiesis:
    return BMHaematopoiesis(port="tcp://*:19010", output_dir=run_issl_dir)


def _make_thymus(run_issl_dir: Path, baseline_override: float | None = None) -> ThymusSelection:
    return ThymusSelection(port="tcp://*:19012", output_dir=run_issl_dir,
                           baseline_import_override=baseline_override)


def _make_transit(run_issl_dir: Path) -> BloodTransit:
    return BloodTransit(port="tcp://*:19011", output_dir=run_issl_dir)


def _make_pln(run_issl_dir: Path) -> PeripheralLN:
    return PeripheralLN(port="tcp://*:19013", output_dir=run_issl_dir)


def _save(record: dict, issl_dir: Path, model_id: str, sim_time_s: float) -> None:
    write_checkpoint(record, issl_dir / model_id, sim_time_s)


# ── Scenario 1: BM baseline (30 days, 6 h steps) ─────────────────────────────

def run_bm1(run_dir: Path) -> None:
    print("  Running BM1 baseline …")
    issl = run_dir / "issl"
    bm = _make_bm(issl)
    t = 0.0
    end = 30 * 86_400
    while t < end:
        bm._step(t, [])
        t += bm.DELTA_T_S
        rec = bm.emit_issl(t)
        _save(rec, issl, "bm_haematopoiesis", t)
    print(f"    BM1 done — export flux at day 30: {rec['export_signals'][0]['flux']:.1f} cells/day")


# ── Scenario 2: THY1 baseline (30 days, 24 h steps) ──────────────────────────

def run_thy1(run_dir: Path) -> None:
    print("  Running THY1 baseline …")
    issl = run_dir / "issl"
    thy = _make_thymus(issl)
    t = 0.0
    end = 30 * 86_400
    while t < end:
        thy._step(t, [])
        t += thy.DELTA_T_S
        rec = thy.emit_issl(t)
        _save(rec, issl, "thymus_selection", t)
    ip = {x["param_id"]: x["value"] for x in rec["internal_parameters"]}
    print(f"    THY1 done — export: {rec['export_signals'][0]['flux']:.1f} agents/cp "
          f"(~{rec['export_signals'][0]['biological_flux_per_day']:,.0f} cells/day scaled)")


# ── Scenario 3: COMP1 — BM → Thymus direct (30 days) ─────────────────────────

def run_comp1(run_dir: Path) -> None:
    print("  Running COMP1 (BM → Thymus direct) …")
    issl = run_dir / "issl"
    bm  = _make_bm(issl)
    thy = _make_thymus(issl, baseline_override=0.0)

    t     = 0.0
    t_thy = 0.0
    end   = 30 * 86_400
    bm_signal: dict | None = None

    while t < end:
        # BM step
        bm._step(t, [])
        t += bm.DELTA_T_S
        bm_rec = bm.emit_issl(t)
        _save(bm_rec, issl, "bm_haematopoiesis", t)
        bm_signal = bm_rec["export_signals"][0]

        # Thymus steps whenever BM has advanced past the next thymus checkpoint
        while t_thy + thy.DELTA_T_S <= t:
            signals = [bm_signal] if bm_signal else []
            thy._step(t_thy, signals)
            t_thy += thy.DELTA_T_S
            thy_rec = thy.emit_issl(t_thy)
            _save(thy_rec, issl, "thymus_selection", t_thy)

    print(f"    COMP1 done")


# ── Scenario 4: COMP2 — BM → Blood Transit → Thymus (30 days) ────────────────

def run_comp2(run_dir: Path) -> None:
    print("  Running COMP2 (BM → Blood Transit → Thymus) …")
    issl = run_dir / "issl"
    bm  = _make_bm(issl)
    bt  = _make_transit(issl)
    thy = _make_thymus(issl, baseline_override=0.0)

    # Pending thymus signals: (deliver_at_s, signal_dict)
    pending: deque[tuple[float, dict]] = deque()

    t     = 0.0
    t_thy = 0.0
    end   = 30 * 86_400

    while t < end:
        bm._step(t, [])
        t += bm.DELTA_T_S
        bm_rec = bm.emit_issl(t)
        _save(bm_rec, issl, "bm_haematopoiesis", t)

        # Run blood transit with BM signal
        bm_sig = bm_rec["export_signals"][0]
        bt._step(t, [bm_sig])
        bt_rec = bt.emit_issl(t)
        _save(bt_rec, issl, "blood_transit", t)

        bt_sig  = bt_rec["export_signals"][0]
        lag_s   = bt_sig.get("lag_s", 0.0) or 0.0
        deliver_at = t + lag_s
        pending.append((deliver_at, bt_sig))

        # Thymus steps
        while t_thy + thy.DELTA_T_S <= t:
            t_thy += thy.DELTA_T_S
            # Collect any signals ready to deliver
            ready = []
            still_pending: deque = deque()
            for (da, sig) in pending:
                if da <= t_thy:
                    ready.append(sig)
                else:
                    still_pending.append((da, sig))
            pending.clear()
            pending.extend(still_pending)

            thy._step(t_thy - thy.DELTA_T_S, ready)
            thy_rec = thy.emit_issl(t_thy)
            _save(thy_rec, issl, "thymus_selection", t_thy)

    print(f"    COMP2 done — transit lag: {lag_s/3600:.1f} h")


# ── Scenario 5: COMP3 — full pipeline (30 days) ───────────────────────────────

def run_comp3(run_dir: Path) -> None:
    print("  Running COMP3 (BM → Transit → Thymus → PLN) …")
    issl = run_dir / "issl"
    bm  = _make_bm(issl)
    bt  = _make_transit(issl)
    thy = _make_thymus(issl, baseline_override=0.0)
    pln = _make_pln(issl)

    # Pending deliveries: (deliver_at_s, signal, target)
    bt_pending:  deque[tuple[float, dict]] = deque()   # → thymus
    thy_pending: deque[tuple[float, dict]] = deque()   # → PLN (2-day lag)

    THY_TO_PLN_LAG_S = 2 * 86_400  # 2-day homing lag

    t     = 0.0
    t_thy = 0.0
    t_pln = 0.0
    end   = 30 * 86_400

    while t < end:
        # BM step (6 h)
        bm._step(t, [])
        t += bm.DELTA_T_S
        bm_rec = bm.emit_issl(t)
        _save(bm_rec, issl, "bm_haematopoiesis", t)

        bm_sig = bm_rec["export_signals"][0]
        bt._step(t, [bm_sig])
        bt_rec = bt.emit_issl(t)
        _save(bt_rec, issl, "blood_transit", t)

        bt_sig  = bt_rec["export_signals"][0]
        lag_s   = bt_sig.get("lag_s", 0.0) or 0.0
        bt_pending.append((t + lag_s, bt_sig))

        # Thymus steps (24 h)
        while t_thy + thy.DELTA_T_S <= t:
            t_thy += thy.DELTA_T_S
            ready_bt: list[dict] = []
            still: deque = deque()
            for (da, sig) in bt_pending:
                (ready_bt if da <= t_thy else still).append(sig) if da <= t_thy else still.append((da, sig))
            bt_pending.clear()
            bt_pending.extend(still)

            thy._step(t_thy - thy.DELTA_T_S, ready_bt)
            thy_rec = thy.emit_issl(t_thy)
            _save(thy_rec, issl, "thymus_selection", t_thy)

            thy_sig = thy_rec["export_signals"][0]
            thy_pending.append((t_thy + THY_TO_PLN_LAG_S, thy_sig))

        # PLN steps (12 h)
        while t_pln + pln.DELTA_T_S <= t:
            t_pln += pln.DELTA_T_S
            ready_thy: list[dict] = []
            still2: deque = deque()
            for item in thy_pending:
                da, sig = item
                if da <= t_pln:
                    ready_thy.append(sig)
                else:
                    still2.append(item)
            thy_pending.clear()
            thy_pending.extend(still2)

            pln._step(t_pln - pln.DELTA_T_S, ready_thy)
            pln_rec = pln.emit_issl(t_pln)
            _save(pln_rec, issl, "peripheral_ln", t_pln)

    print(f"    COMP3 done")


# ── Figure generation ─────────────────────────────────────────────────────────

def generate_figures() -> None:
    print("\nGenerating figures …")
    import subprocess
    repo = Path(__file__).parent
    venv_python = repo / ".venv" / "bin" / "python"
    python = str(venv_python) if venv_python.exists() else sys.executable
    result = subprocess.run(
        [python, "results/plot_results.py"],
        cwd=str(repo),
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  plot_results.py stderr:", result.stderr[:500])
    print(result.stdout.strip())


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figs", action="store_true", help="figures only")
    args = ap.parse_args()

    if not args.figs:
        print("Running all scenarios …")
        scenarios = [
            ("BM1_baseline",    run_bm1),
            ("THY1_baseline",   run_thy1),
            ("COMP1_direct",    run_comp1),
            ("COMP2_transfer",  run_comp2),
            ("COMP3_full_graph", run_comp3),
        ]
        for name, fn in scenarios:
            run_dir = RESULTS / name / name
            run_dir.mkdir(parents=True, exist_ok=True)
            fn(run_dir)
        print("\nAll simulations complete.")

    generate_figures()


if __name__ == "__main__":
    main()
