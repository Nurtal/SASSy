#!/usr/bin/env python3
"""
Model 4 — Peripheral Lymph Node Dynamics (ODE)

State vector: [CD4, CD8]  (naïve T-cell pool sizes)

ODEs (Borghans-De Boer homeostasis model):
  dCD4/dt = import_CD4 + rho_4 * (S4 - CD4) - (d_act + d_death) * CD4
  dCD8/dt = import_CD8 + rho_8 * (S8 - CD8) - (d_act + d_death) * CD8

Import signals are received as piecewise constant terms updated each step.

Usage:
  python models/peripheral_ln/model.py --port tcp://*:5013 --output-dir logs/run/issl
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

_ID_CD4 = "CL:0000624"  # CD4+ naive T cell
_ID_CD8 = "CL:0000625"  # CD8+ naive T cell


def _load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


class PeripheralLN(ModelBase):
    """Peripheral lymph node ODE model."""

    MODEL_ID = "peripheral_ln"
    MODEL_VERSION = "1"
    DELTA_T_S = 43_200.0  # 12 hours

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
        self._state = np.array([ic["CD4"], ic["CD8"]], dtype=float)

        # Import rates (cells/day) — updated by _step when signals arrive
        self._import_cd4: float = 0.0
        self._import_cd8: float = 0.0
        self._last_import_rate: float = 0.0

    # ------------------------------------------------------------------
    # ODE system

    def _ode(self, t: float, y: np.ndarray) -> np.ndarray:
        CD4, CD8 = y
        p = self._p
        dCD4 = (self._import_cd4
                + p["rho_4"] * (p["S4"] - CD4)
                - (p["d_act"] + p["d_death"]) * CD4)
        dCD8 = (self._import_cd8
                + p["rho_8"] * (p["S8"] - CD8)
                - (p["d_act"] + p["d_death"]) * CD8)
        return np.array([dCD4, dCD8])

    # ------------------------------------------------------------------
    # Abstract interface

    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        """Advance the ODE by one DELTA_T_S, incorporating any import signals."""
        # Update import rates from incoming signals
        for sig in signals:
            if "flux" in sig:
                flux_per_checkpoint = float(sig["flux"])
                # flux is cells·checkpoint^-1 from thymus; convert to cells/day
                flux_per_day = flux_per_checkpoint  # thymus emits per 24-h checkpoint
                self._last_import_rate = flux_per_day
                self._import_cd4 = flux_per_day * self._p["cd4_fraction"]
                self._import_cd8 = flux_per_day * (1.0 - self._p["cd4_fraction"])
                break

        t_start_d = sim_time_s / 86_400
        t_end_d = (sim_time_s + self.DELTA_T_S) / 86_400

        sol = solve_ivp(
            self._ode,
            [t_start_d, t_end_d],
            self._state,
            method="RK45",
            rtol=1e-6,
            atol=1e-8,
        )
        if not sol.success:
            logger.warning("ODE integration warning at t=%s: %s", sim_time_s, sol.message)

        self._state = np.maximum(sol.y[:, -1], 0.0)

    def emit_issl(self, sim_time_s: float) -> dict:
        CD4, CD8 = self._state
        total = CD4 + CD8
        ratio = CD4 / CD8 if CD8 > 0 else float("inf")

        ci_cd4 = self._state_ci95(CD4, "rho_4", "S4")
        ci_cd8 = self._state_ci95(CD8, "rho_8", "S8")
        homeo_rate = self._p["rho_4"] * (self._p["S4"] - CD4)

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
                    "entity_id":    _ID_CD4,
                    "label":        "Naïve CD4+ T cells",
                    "count":        round(CD4, 1),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD4+", "CD8-", "CD44lo", "CD62Lhi"],
                    "ci_95":        ci_cd4,
                },
                {
                    "entity_class": "obo",
                    "entity_id":    _ID_CD8,
                    "label":        "Naïve CD8+ T cells",
                    "count":        round(CD8, 1),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD4-", "CD8+", "CD44lo", "CD62Lhi"],
                    "ci_95":        ci_cd8,
                },
            ],
            "discrete_events": [],
            # No export signals in COMP-3 baseline configuration
            "export_signals": [],
            "internal_parameters": [
                {
                    "param_id": k,
                    "value":    round(v, 6),
                    "unit":     self._p_meta[k]["unit"],
                    "ci_95":    self._p_ci[k],
                    "source":   self._p_meta[k].get("source", ""),
                }
                for k, v in self._p.items()
            ],
            "watchdog": self._make_watchdog(
                sim_time_s,
                ood_flag=bool(ratio < 0.5 or ratio > 10.0 or CD4 < 0 or CD8 < 0),
                divergence_score=round(abs(ratio - 2.0) / 2.0, 4),
            ),
        }

    # ------------------------------------------------------------------
    # CI-95 helper

    def _state_ci95(self, value: float, *param_keys: str) -> list[float]:
        variance = 0.0
        for k in param_keys:
            ci = self._p_ci[k]
            sigma = (ci[1] - ci[0]) / (2 * 1.96)
            variance += sigma ** 2
        sigma_total = math.sqrt(variance) * abs(value) * 0.05  # 5% relative sensitivity
        return [round(value - 1.96 * sigma_total, 2), round(value + 1.96 * sigma_total, 2)]


# ------------------------------------------------------------------
# Entry point

def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(name)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="OISA — Peripheral Lymph Node model")
    parser.add_argument("--port", default="tcp://*:5013")
    parser.add_argument("--output-dir", default="logs/run/issl")
    args = parser.parse_args()

    model = PeripheralLN(port=args.port, output_dir=Path(args.output_dir))
    model.run()


if __name__ == "__main__":
    main()
