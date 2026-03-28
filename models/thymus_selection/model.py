#!/usr/bin/env python3
"""
Model 3 — Thymic T-cell Selection (ABM)

Agents: individual thymocytes with a TCR affinity score.
Zones:  cortex (positive selection) + medulla (negative selection).

Per checkpoint (24 h):
  - 12 independent realisations are run
  - Each realisation advances by substep_h-hour steps
  - CI-95 computed empirically across realisations

Selection rules:
  - DN → DP: probabilistic transition in cortex
  - Positive selection: affinity in [theta_low, theta_high]
  - Negative selection: affinity > theta_negative (medulla)
  - Neglect death: DP that never encounters a valid pMHC → dies

Usage:
  python models/thymus_selection/model.py --port tcp://*:5012 --output-dir logs/run/issl
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from models._base.model_base import ModelBase

logger = logging.getLogger(__name__)

_PARAMS_PATH = Path(__file__).parent / "parameters.yaml"

_ID_DN1  = "CL:0002420"  # early T-lineage progenitor
_ID_DP   = "CL:0000893"  # double-positive thymocyte
_ID_CD4SP= "CL:0000624"  # CD4+ single-positive
_ID_CD8SP= "CL:0000625"  # CD8+ single-positive
_ID_NAIVE= "CL:0000898"  # naive T cell (exported)

# OBO IDs for discrete events
_EV_POS_SEL  = "GO:0045058"  # T cell positive selection
_EV_NEG_SEL  = "GO:0045059"  # T cell negative selection
_EV_NEGLECT  = "GO:0070227"  # lymphocyte clonal deletion (neglect reuse)


def _load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


# ---------------------------------------------------------------------------
# Agent

@dataclass
class Thymocyte:
    stage: str          # DN1|DN2|DN3|DN4|DP|CD4SP|CD8SP
    tcr_affinity: float
    zone: str           # cortex|medulla
    age_steps: int = 0
    fate: str = "alive" # alive|pos_selected|neg_deleted|neglect_dead|exported
    lineage: str = ""   # CD4|CD8 (set at DP→SP)
    medullary_dwell: int = 0  # steps spent in medulla after positive selection
    dp_steps: int = 0   # steps spent at DP stage (for max_dp_age_steps check)


# ---------------------------------------------------------------------------
# Single realisation

class ThymicRealisation:
    """One ABM realisation."""

    def __init__(self, p: dict, tcr_cfg: dict, rng: np.random.Generator) -> None:
        self._p = p
        self._tcr_mu = tcr_cfg["mu"]
        self._tcr_sigma = tcr_cfg["sigma"]
        self._rng = rng
        self.agents: list[Thymocyte] = []
        self.exported_this_checkpoint: int = 0
        self.pos_sel_count: int = 0
        self.neg_sel_count: int = 0
        self.neglect_count: int = 0

    def add_progenitors(self, n: int) -> None:
        """Add *n* new DN1 agents with sampled TCR affinities."""
        if n <= 0:
            return
        affinities = self._rng.lognormal(
            mean=self._tcr_mu, sigma=self._tcr_sigma, size=n
        )
        for aff in affinities:
            self.agents.append(Thymocyte(stage="DN1", tcr_affinity=float(aff),
                                          zone="cortex"))

    def step(self) -> None:
        """Advance all agents by one sub-step."""
        p = self._p
        survivors = []
        for agent in self.agents:
            if agent.fate != "alive":
                continue
            agent.age_steps += 1

            # --- DN stages → DP transition in cortex
            if agent.stage in ("DN1", "DN2", "DN3", "DN4"):
                # Simple probabilistic progression through DN stages
                if self._rng.random() < p["p_transition_DP"]:
                    agent.stage = "DP"
                survivors.append(agent)
                continue

            # --- DP stage: positive selection in cortex
            if agent.stage == "DP":
                if agent.zone == "cortex":
                    agent.dp_steps += 1
                    # Age-based neglect death: DP cells that exceed their
                    # biological lifespan without receiving a selection signal die.
                    if agent.dp_steps >= int(p["max_dp_age_steps"]):
                        agent.fate = "neglect_dead"
                        self.neglect_count += 1
                        # Not added to survivors → removed from simulation
                    elif self._rng.random() < p["p_encounter_cortex"]:
                        aff = agent.tcr_affinity
                        if aff < p["theta_low"]:
                            # Neglect death on pMHC encounter (affinity too low)
                            agent.fate = "neglect_dead"
                            self.neglect_count += 1
                        elif aff <= p["theta_high"]:
                            # Positive selection → become SP.
                            # fate stays "alive" so the agent is not skipped by the
                            # `if agent.fate != "alive": continue` guard next sub-step.
                            agent.fate = "alive"
                            self.pos_sel_count += 1
                            lineage = "CD4" if self._rng.random() < p["cd4_fraction"] else "CD8"
                            agent.lineage = lineage
                            agent.stage = "CD4SP" if lineage == "CD4" else "CD8SP"
                            agent.zone = "cortex"  # will migrate to medulla next steps
                        else:
                            # High affinity DP → negative selection at DP stage
                            agent.fate = "neg_deleted"
                            self.neg_sel_count += 1
                        # Only keep agent if still alive (fixes ghost accumulation:
                        # previously survivors.append was unconditional for DP)
                        if agent.fate == "alive":
                            survivors.append(agent)
                    else:
                        # No encounter this step — keep agent alive for next step
                        survivors.append(agent)
                continue

            # --- SP stage: migrate to medulla, then negative selection
            if agent.stage in ("CD4SP", "CD8SP"):
                # Migrate to medulla if still in cortex
                if agent.zone == "cortex":
                    if self._rng.random() < p["p_migrate_medulla"]:
                        agent.zone = "medulla"
                        agent.fate = "alive"  # reset from pos_selected label
                    survivors.append(agent)
                    continue

                # In medulla: negative selection encounters
                agent.medullary_dwell += 1
                if self._rng.random() < p["p_encounter_medulla"]:
                    if agent.tcr_affinity > p["theta_negative"]:
                        agent.fate = "neg_deleted"
                        self.neg_sel_count += 1
                        # Don't keep in survivors
                        continue

                # Export after dwell time
                if agent.medullary_dwell >= int(p["medullary_dwell_steps"]):
                    agent.fate = "exported"
                    self.exported_this_checkpoint += 1
                    # Don't keep in survivors
                    continue

                survivors.append(agent)

        self.agents = survivors

    def reset_checkpoint_counters(self) -> None:
        self.exported_this_checkpoint = 0
        self.pos_sel_count = 0
        self.neg_sel_count = 0
        self.neglect_count = 0

    # Count helpers
    def count_stage(self, stage: str) -> int:
        return sum(1 for a in self.agents if a.stage == stage)

    def count_zone(self, zone: str) -> int:
        return sum(1 for a in self.agents if a.zone == zone and a.fate == "alive")


# ---------------------------------------------------------------------------
# Model

class ThymusSelection(ModelBase):
    """Thymic T-cell selection ABM."""

    MODEL_ID = "thymus_selection"
    MODEL_VERSION = "3"
    DELTA_T_S = 86_400.0  # 24 h

    def __init__(self, port: str, output_dir: Path,
                 baseline_import_override: float | None = None) -> None:
        super().__init__(
            model_id=self.MODEL_ID,
            model_version=self.MODEL_VERSION,
            formalism="ABM",
            delta_t_s=self.DELTA_T_S,
            port=port,
            output_dir=output_dir,
        )
        cfg = _load_params()
        self._p = {k: v["value"] for k, v in cfg["parameters"].items()}
        self._p_ci = {k: v["ci_95"] for k, v in cfg["parameters"].items()}
        self._p_meta = cfg["parameters"]
        self._tcr_cfg = cfg["tcr_affinity"]
        self._scaling_factor: float = float(self._p.get("scaling_factor", 1.0))
        self._n_realisations: int = cfg["abm"]["n_realisations"]
        self._substep_h: int = cfg["abm"]["substep_h"]
        self._steps_per_checkpoint = int(24 / self._substep_h)

        # Allow caller to override baseline_import (e.g. pass 0.0 in coupled
        # runs to prevent artificial pool pre-loading during upstream lag periods).
        if baseline_import_override is not None:
            self._p["baseline_import"] = baseline_import_override

        # One RNG per realisation for reproducibility
        base_rng = np.random.default_rng(seed=42)
        self._realisations: list[ThymicRealisation] = [
            ThymicRealisation(
                p=self._p,
                tcr_cfg=self._tcr_cfg,
                rng=np.random.default_rng(base_rng.integers(0, 2**32)),
            )
            for _ in range(self._n_realisations)
        ]

        # Results of the last checkpoint (means + ci_95 across realisations)
        self._last_sim_time_s: float = 0.0
        self._last_results: dict = {}

    # ------------------------------------------------------------------
    # Abstract interface

    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        """Advance all realisations by one checkpoint (24 h = substep_h-h sub-steps)."""
        self._last_sim_time_s = sim_time_s + self.DELTA_T_S

        # Parse incoming progenitor signal.
        # If no upstream edge delivers a signal (e.g. THY1 standalone), fall back
        # to the baseline_import parameter (steady-state BM export estimate).
        signal_received = False
        n_import_mean = 0.0
        n_import_sigma = 0.0
        for sig in signals:
            if "flux" in sig:
                flux = float(sig["flux"])
                n_import_mean = max(flux, 0.0)
                if "ci_95" in sig:
                    ci = sig["ci_95"]
                    n_import_sigma = (ci[1] - ci[0]) / (2 * 1.96)
                signal_received = True
                break
        if not signal_received:
            baseline = self._p.get("baseline_import", 0.0)
            if baseline > 0.0:
                n_import_mean = baseline
                ci = self._p_ci.get("baseline_import", [baseline, baseline])
                n_import_sigma = (ci[1] - ci[0]) / (2 * 1.96)

        # Per-realisation results
        exported_counts = []
        pos_sel_counts = []
        neg_sel_counts = []
        neglect_counts = []
        dn_counts = []
        dp_counts = []
        cd4sp_counts = []
        cd8sp_counts = []

        for real in self._realisations:
            real.reset_checkpoint_counters()

            # Sample new arrivals for this realisation
            n_new = int(round(max(
                real._rng.normal(n_import_mean, max(n_import_sigma, 1.0)),
                0.0
            )))
            real.add_progenitors(n_new)

            # Run sub-steps
            for _ in range(self._steps_per_checkpoint):
                real.step()

            exported_counts.append(real.exported_this_checkpoint)
            pos_sel_counts.append(real.pos_sel_count)
            neg_sel_counts.append(real.neg_sel_count)
            neglect_counts.append(real.neglect_count)
            dn_counts.append(sum(real.count_stage(s) for s in ("DN1","DN2","DN3","DN4")))
            dp_counts.append(real.count_stage("DP"))
            cd4sp_counts.append(real.count_stage("CD4SP"))
            cd8sp_counts.append(real.count_stage("CD8SP"))

        def _stats(vals: list[int | float]) -> tuple[float, list[float]]:
            arr = np.array(vals, dtype=float)
            mean = float(np.mean(arr))
            lo, hi = float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))
            return mean, [round(lo, 2), round(hi, 2)]

        self._last_results = {
            "exported":  _stats(exported_counts),
            "pos_sel":   _stats(pos_sel_counts),
            "neg_sel":   _stats(neg_sel_counts),
            "neglect":   _stats(neglect_counts),
            "dn":        _stats(dn_counts),
            "dp":        _stats(dp_counts),
            "cd4sp":     _stats(cd4sp_counts),
            "cd8sp":     _stats(cd8sp_counts),
        }

    def emit_issl(self, sim_time_s: float) -> dict:
        r = self._last_results
        if not r:
            # No step has run yet — return zero-state
            zero = (0.0, [0.0, 0.0])
            r = {k: zero for k in ("exported","pos_sel","neg_sel","neglect",
                                    "dn","dp","cd4sp","cd8sp")}

        n_real = self._n_realisations
        exported_mean, exported_ci = r["exported"]
        sf = self._scaling_factor

        def _scale(val: float) -> float:
            return round(val * sf, 1)

        def _scale_ci(ci: list[float]) -> list[float]:
            return [round(ci[0] * sf, 1), round(ci[1] * sf, 1)]

        return {
            "envelope": {
                "model_id":      self.MODEL_ID,
                "model_version": self.MODEL_VERSION,
                "sim_time_s":    sim_time_s,
                "schema_uri":    "schemas/issl_v1.schema.json",
                "formalism":     "ABM",
            },
            # counts are scaled (ABM agent count × scaling_factor) to give
            # biologically realistic cell numbers.  Raw agent counts are
            # preserved in internal_parameters (raw_*_agents fields).
            "continuous_state": [
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_DN1,
                    "label":        "Double-negative thymocytes (DN1-DN4)",
                    "count":        _scale(r["dn"][0]),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD44+", "CD25+/-", "CD117+"],
                    "ci_95":        _scale_ci(r["dn"][1]),
                },
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_DP,
                    "label":        "Double-positive thymocytes",
                    "count":        _scale(r["dp"][0]),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD4+", "CD8+", "TCRαβ+"],
                    "ci_95":        _scale_ci(r["dp"][1]),
                },
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_CD4SP,
                    "label":        "CD4+ single-positive thymocytes",
                    "count":        _scale(r["cd4sp"][0]),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD4+", "CD8-", "TCRαβhi"],
                    "ci_95":        _scale_ci(r["cd4sp"][1]),
                },
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_CD8SP,
                    "label":        "CD8+ single-positive thymocytes",
                    "count":        _scale(r["cd8sp"][0]),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD4-", "CD8+", "TCRαβhi"],
                    "ci_95":        _scale_ci(r["cd8sp"][1]),
                },
            ],
            "discrete_events": [
                {
                    "event_type":     "positive_selection",
                    "entity_id":      _EV_POS_SEL,
                    "count":          int(round(r["pos_sel"][0])),
                    "sim_time_s":     sim_time_s,
                    "n_realisations": n_real,
                    "variance":       self._empirical_variance(r["pos_sel"]),
                },
                {
                    "event_type":     "negative_selection",
                    "entity_id":      _EV_NEG_SEL,
                    "count":          int(round(r["neg_sel"][0])),
                    "sim_time_s":     sim_time_s,
                    "n_realisations": n_real,
                    "variance":       self._empirical_variance(r["neg_sel"]),
                },
                {
                    "event_type":     "neglect_death",
                    "entity_id":      _EV_NEGLECT,
                    "count":          int(round(r["neglect"][0])),
                    "sim_time_s":     sim_time_s,
                    "n_realisations": n_real,
                    "variance":       self._empirical_variance(r["neglect"]),
                },
            ],
            "export_signals": [
                {
                    # flux is at RAW ABM scale (agents·checkpoint⁻¹) for
                    # orchestrator / PLN coupling.
                    # Biological estimate = flux × scaling_factor.
                    "signal_id":              "thymus_selection.naive_T_export",
                    "entity_id":              _ID_NAIVE,
                    "flux":                   round(exported_mean, 2),
                    "unit":                   "agents·checkpoint^-1",
                    "lag_s":                  None,
                    "ci_95":                  exported_ci,
                    "scaling_factor":         sf,
                    "biological_flux_per_day": round(exported_mean * sf, 0),
                    "biological_unit":        "cells·day^-1",
                }
            ],
            "internal_parameters": [
                {
                    "param_id": k,
                    "value":    v,
                    "unit":     self._p_meta[k]["unit"],
                    "ci_95":    self._p_ci[k],
                    "source":   self._p_meta[k].get("source", ""),
                }
                for k, v in self._p.items()
            ] + [
                # Raw ABM agent counts — multiply by scaling_factor for biological estimate
                {
                    "param_id": "raw_dn_agents",
                    "value":    round(r["dn"][0], 1),
                    "unit":     "agents",
                    "ci_95":    r["dn"][1],
                    "source":   "ABM realisation mean",
                },
                {
                    "param_id": "raw_dp_agents",
                    "value":    round(r["dp"][0], 1),
                    "unit":     "agents",
                    "ci_95":    r["dp"][1],
                    "source":   "ABM realisation mean",
                },
                {
                    "param_id": "raw_cd4sp_agents",
                    "value":    round(r["cd4sp"][0], 1),
                    "unit":     "agents",
                    "ci_95":    r["cd4sp"][1],
                    "source":   "ABM realisation mean",
                },
                {
                    "param_id": "raw_cd8sp_agents",
                    "value":    round(r["cd8sp"][0], 1),
                    "unit":     "agents",
                    "ci_95":    r["cd8sp"][1],
                    "source":   "ABM realisation mean",
                },
                {
                    "param_id": "raw_exported_agents",
                    "value":    round(exported_mean, 2),
                    "unit":     "agents·checkpoint^-1",
                    "ci_95":    exported_ci,
                    "source":   "ABM realisation mean",
                },
            ],
            "watchdog": self._make_watchdog(
                sim_time_s,
                ood_flag=bool(exported_mean < 0),
                divergence_score=None,
            ),
        }

    # ------------------------------------------------------------------
    # Helpers

    @staticmethod
    def _empirical_variance(stats: tuple[float, list[float]]) -> float:
        """Approximate variance from empirical CI-95."""
        mean, (lo, hi) = stats
        sigma = (hi - lo) / (2 * 1.96)
        return round(sigma ** 2, 4)


# ------------------------------------------------------------------
# Entry point

def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(name)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="OISA — Thymus Selection ABM")
    parser.add_argument("--port", default="tcp://*:5012")
    parser.add_argument("--output-dir", default="logs/run/issl")
    parser.add_argument(
        "--baseline-import", type=float, default=None,
        help="Override baseline progenitor import (cells/day). "
             "Pass 0 in coupled runs to prevent artificial pre-loading during upstream lags. "
             "Defaults to the baseline_import value in parameters.yaml.",
    )
    args = parser.parse_args()

    model = ThymusSelection(
        port=args.port,
        output_dir=Path(args.output_dir),
        baseline_import_override=args.baseline_import,
    )
    model.run()


if __name__ == "__main__":
    main()
