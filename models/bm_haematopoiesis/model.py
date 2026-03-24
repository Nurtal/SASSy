#!/usr/bin/env python3
"""
Model 1 — Bone Marrow Haematopoiesis (ODE)

State vector: [HSC, CLP, DN1]

ODEs:
  dHSC/dt = r_self * HSC * (1 - HSC/K_niche) - d_diff * HSC - d_apop * HSC
  dCLP/dt = d_diff * HSC - (alpha_T + alpha_myeloid) * CLP - d_CLP * CLP
  dDN1/dt = alpha_T * CLP - export_rate * DN1

Export flux (cells/day): export_rate * DN1
CI-95: propagated analytically from parameter posterior variances.

Usage:
  python models/bm_haematopoiesis/model.py \
      --port tcp://*:5010 \
      --output-dir logs/run/issl
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from models._base.model_base import ModelBase

logger = logging.getLogger(__name__)

_PARAMS_PATH = Path(__file__).parent / "parameters.yaml"

# OBO entity IDs
_ID_HSC = "CL:0000037"
_ID_CLP = "CL:0000051"
_ID_DN1 = "CL:0002420"  # early T-lineage progenitor / DN1/ETP


def _load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


class BMHaematopoiesis(ModelBase):
    """Bone marrow ODE model."""

    MODEL_ID = "bm_haematopoiesis"
    MODEL_VERSION = "3"
    DELTA_T_S = 21_600.0  # 6 hours in seconds

    def __init__(self, port: str, output_dir: Path) -> None:
        super().__init__(
            model_id=self.MODEL_ID,
            model_version=self.MODEL_VERSION,
            formalism="ODE",
            delta_t_s=self.DELTA_T_S,
            port=port,
            output_dir=output_dir,
        )
        cfg = _load_params()
        self._p = {k: v["value"] for k, v in cfg["parameters"].items()}
        self._p_ci = {k: v["ci_95"] for k, v in cfg["parameters"].items()}
        self._p_meta = cfg["parameters"]

        ic = cfg["initial_conditions"]
        # State: [HSC, CLP, DN1]
        self._state = np.array([ic["HSC"], ic["CLP"], ic["DN1"]], dtype=float)
        self._sim_time_s: float = 0.0

    # ------------------------------------------------------------------
    # ODE system

    def _ode(self, t: float, y: np.ndarray) -> np.ndarray:
        HSC, CLP, DN1 = y
        p = self._p
        dHSC = p["r_self"] * HSC * (1 - HSC / p["K_niche"]) - p["d_diff"] * HSC - p["d_apop"] * HSC
        dCLP = p["d_diff"] * HSC - (p["alpha_T"] + p["alpha_myeloid"]) * CLP - p["d_CLP"] * CLP
        dDN1 = p["alpha_T"] * CLP - p["export_rate"] * DN1
        return np.array([dHSC, dCLP, dDN1])

    # ------------------------------------------------------------------
    # Abstract interface

    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        """Advance the ODE by one delta_t_s starting from *sim_time_s*."""
        # Niche cytokine signal (optional — currently unused in baseline)
        _ = signals

        t_start = sim_time_s
        t_end = sim_time_s + self.DELTA_T_S
        # Convert from seconds to days for the ODE (parameters are in day^-1)
        t_start_d = t_start / 86_400
        t_end_d = t_end / 86_400

        sol = solve_ivp(
            self._ode,
            [t_start_d, t_end_d],
            self._state,
            method="RK45",
            dense_output=False,
            rtol=1e-6,
            atol=1e-8,
        )
        if not sol.success:
            logger.warning("ODE integration warning at t=%s: %s", sim_time_s, sol.message)

        self._state = np.maximum(sol.y[:, -1], 0.0)  # clamp negatives from numerics
        self._sim_time_s = t_end

    def emit_issl(self, sim_time_s: float) -> dict:
        HSC, CLP, DN1 = self._state
        p = self._p
        export_flux = p["export_rate"] * DN1  # cells/day

        # Analytical CI-95 on export flux via first-order sensitivity
        # Var(export_rate * DN1) ≈ DN1^2 * Var(export_rate) + export_rate^2 * Var(DN1)
        # Var(DN1) approximated from parameter CI via sensitivity coefficients
        ci_export = self._flux_ci95(DN1, export_flux)
        ci_hsc = self._state_ci95(HSC, "r_self", "K_niche")
        ci_clp = self._state_ci95(CLP, "alpha_T", "alpha_myeloid")
        ci_dn1 = self._state_ci95(DN1, "alpha_T", "export_rate")

        return {
            "envelope": {
                "model_id":      self.MODEL_ID,
                "model_version": self.MODEL_VERSION,
                "sim_time_s":    sim_time_s,
                "schema_uri":    "schemas/issl_v1.schema.json",
                "formalism":     "ODE",
            },
            "continuous_state": [
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_HSC,
                    "label":        "Haematopoietic stem cell",
                    "count":        round(HSC, 2),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD34+", "CD117+", "CD150+", "CD48-"],
                    "ci_95":        ci_hsc,
                },
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_CLP,
                    "label":        "Common lymphoid progenitor",
                    "count":        round(CLP, 2),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD34+", "CD127+", "CD117lo"],
                    "ci_95":        ci_clp,
                },
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_DN1,
                    "label":        "Early T-lineage progenitor (DN1/ETP)",
                    "count":        round(DN1, 2),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD34+", "CD117+", "CD44+", "CD25-"],
                    "ci_95":        ci_dn1,
                },
            ],
            "discrete_events": [],
            "export_signals": [
                {
                    "signal_id": "bm_haematopoiesis.progenitor_export",
                    "entity_id": _ID_DN1,
                    "flux":      round(export_flux, 4),
                    "unit":      "cells·day^-1",
                    "lag_s":     None,
                    "ci_95":     ci_export,
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
                for k, v in p.items()
            ],
            "watchdog": self._make_watchdog(
                sim_time_s,
                ood_flag=self._is_ood(HSC, CLP, DN1),
                divergence_score=self._divergence(HSC, CLP, DN1),
            ),
        }

    # ------------------------------------------------------------------
    # CI-95 helpers (first-order sensitivity / diagonal covariance)

    def _state_ci95(self, value: float, *param_keys: str) -> list[float]:
        """Approximate CI-95 for a state variable based on parameter posteriors."""
        variance = 0.0
        for k in param_keys:
            ci = self._p_ci[k]
            # Treat CI as ±1.96σ → σ = (upper - lower) / (2 * 1.96)
            sigma = (ci[1] - ci[0]) / (2 * 1.96)
            variance += sigma ** 2
        sigma_total = math.sqrt(variance) * value
        return [round(value - 1.96 * sigma_total, 2), round(value + 1.96 * sigma_total, 2)]

    def _flux_ci95(self, dn1: float, flux: float) -> list[float]:
        p = self._p
        ci_er = self._p_ci["export_rate"]
        sigma_er = (ci_er[1] - ci_er[0]) / (2 * 1.96)
        sigma_flux = math.sqrt((dn1 * sigma_er) ** 2)
        return [round(flux - 1.96 * sigma_flux, 4), round(flux + 1.96 * sigma_flux, 4)]

    def _is_ood(self, HSC: float, CLP: float, DN1: float) -> bool:
        return bool(HSC < 500 or HSC > 50_000 or CLP < 0 or DN1 < 0)

    def _divergence(self, HSC: float, CLP: float, DN1: float) -> float:
        """Scalar divergence score: relative deviation of HSC from niche set point."""
        return abs(HSC - self._p["K_niche"]) / self._p["K_niche"]


# ------------------------------------------------------------------
# Entry point

def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(name)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="OISA — BM Haematopoiesis model")
    parser.add_argument("--port", default="tcp://*:5010",
                        help="ZMQ REP bind address (default: tcp://*:5010)")
    parser.add_argument("--output-dir", default="logs/run/issl",
                        help="Base directory for ISSL checkpoint files")
    args = parser.parse_args()

    model = BMHaematopoiesis(port=args.port, output_dir=Path(args.output_dir))
    model.run()


if __name__ == "__main__":
    main()
