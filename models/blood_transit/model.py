#!/usr/bin/env python3
"""
Model 2 — Blood Transit Kinetics (ODE, transfer model)

This is a *stateless* transfer model: it is invoked on-demand by the
orchestrator when a BM progenitor export signal arrives on the
BM→Thymus edge.  It does not maintain persistent state across global
checkpoints.

Each invocation:
  1. Receives N0 = incoming signal flux as initial condition
  2. Integrates dN/dt = -(k_homing + k_death) * N until 95 % of N0 is gone
  3. Computes cells_delivered = ∫ k_homing * N(t) dt
  4. Reports cells_delivered and lag_s in export_signals

The lag_s field is read by the orchestrator to schedule delivery of
progenitors to the thymus model.

Usage (as transfer model, invoked by orchestrator — not run continuously):
  python models/blood_transit/model.py --port tcp://*:5011 --output-dir logs/run/issl
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

_ID_DN1 = "CL:0002420"


def _load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


class BloodTransit(ModelBase):
    """Blood transit ODE — transfer model."""

    MODEL_ID = "blood_transit"
    MODEL_VERSION = "1"

    def __init__(self, port: str, output_dir: Path) -> None:
        # delta_t_s=None because this model is invoked on-demand, not on a clock
        super().__init__(
            model_id=self.MODEL_ID,
            model_version=self.MODEL_VERSION,
            formalism="ODE",
            delta_t_s=None,
            port=port,
            output_dir=output_dir,
        )
        cfg = _load_params()
        self._p = {k: v["value"] for k, v in cfg["parameters"].items()}
        self._p_ci = {k: v["ci_95"] for k, v in cfg["parameters"].items()}
        self._p_meta = cfg["parameters"]

        # Results of the last transit computation
        self._last_cells_delivered: float = 0.0
        self._last_lag_s: float = 0.0
        self._last_viability: float = 0.0
        self._last_n0: float = 0.0
        self._last_sim_time_s: float = 0.0

    # ------------------------------------------------------------------
    # ODE system

    def _ode(self, t: float, y: np.ndarray) -> np.ndarray:
        """First-order transit: dN/dt = -(k_homing + k_death) * N"""
        N = y[0]
        loss_rate = self._p["k_homing"] + self._p["k_death"]
        return np.array([-loss_rate * N])

    # ------------------------------------------------------------------
    # Abstract interface

    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        """Run transit ODE from the incoming signal's flux as N0.

        Unlike clock-driven models, this runs the full transit simulation
        immediately and stores results for emit_issl().
        """
        if not signals:
            logger.warning("BloodTransit._step called with no incoming signal — skipping.")
            self._last_n0 = 0.0
            self._last_cells_delivered = 0.0
            self._last_lag_s = 0.0
            self._last_viability = 0.0
            self._last_sim_time_s = sim_time_s
            return

        # Take the first signal (BM progenitor export)
        signal = signals[0]
        n0 = float(signal.get("flux", 0.0))
        self._last_n0 = n0
        self._last_sim_time_s = sim_time_s

        if n0 <= 0:
            self._last_cells_delivered = 0.0
            self._last_lag_s = 0.0
            self._last_viability = 0.0
            return

        k_homing = self._p["k_homing"]
        k_death = self._p["k_death"]
        k_total = k_homing + k_death
        stop_frac = self._p["stop_fraction"]

        # Analytical solution: N(t) = n0 * exp(-k_total * t_days)
        # Cumulative homed: H(t) = (k_homing / k_total) * n0 * (1 - exp(-k_total * t))
        # Cumulative lost:  L(t) = (k_death / k_total) * n0 * (1 - exp(-k_total * t))
        # H(t) + L(t) = n0 * (1 - exp(-k_total * t))
        # Stop when H + L >= stop_frac * n0:
        #   exp(-k_total * t) = 1 - stop_frac
        #   t_stop = -ln(1 - stop_frac) / k_total  (in days)

        t_stop_days = -math.log(1.0 - stop_frac) / k_total

        # Integrate numerically to get cumulative homing flux
        def augmented_ode(t, y):
            N, cum_homed = y
            dN = -k_total * N
            d_homed = k_homing * N
            return [dN, d_homed]

        sol = solve_ivp(
            augmented_ode,
            [0.0, t_stop_days],
            [n0, 0.0],
            method="RK45",
            rtol=1e-8,
            atol=1e-10,
        )

        cells_delivered = float(sol.y[1, -1])
        lag_days = sol.t[-1]
        lag_s = lag_days * 86_400

        self._last_cells_delivered = max(cells_delivered, 0.0)
        self._last_lag_s = lag_s
        self._last_viability = cells_delivered / n0 if n0 > 0 else 0.0

    def emit_issl(self, sim_time_s: float) -> dict:
        n0 = self._last_n0

        # CI-95 on delivered cells from parameter uncertainty
        ci_delivered = self._delivered_ci95(n0, self._last_cells_delivered)
        ci_lag = self._lag_ci95(self._last_lag_s)

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
                    "entity_id":    _ID_DN1,
                    "label":        "DN1 progenitors in transit (blood)",
                    "count":        round(self._last_cells_delivered, 4),
                    "unit":         "cells",
                    "fitness":      None,
                    "surface_markers": ["CD34+", "CD117+", "CD44+"],
                    "ci_95":        ci_delivered,
                }
            ],
            "discrete_events": [],
            "export_signals": [
                {
                    "signal_id": "blood_transit.delivered",
                    "entity_id": _ID_DN1,
                    "flux":      round(self._last_cells_delivered, 4),
                    "unit":      "cells",
                    # lag_s is the key field read by the orchestrator's TransferDispatcher
                    "lag_s":     round(self._last_lag_s, 2),
                    "ci_95":     ci_delivered,
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
            ],
            "watchdog": self._make_watchdog(
                sim_time_s,
                ood_flag=bool(self._last_viability < 0.1 and n0 > 0),
                divergence_score=round(1.0 - self._last_viability, 4) if n0 > 0 else None,
            ),
        }

    # ------------------------------------------------------------------
    # CI-95 helpers

    def _delivered_ci95(self, n0: float, delivered: float) -> list[float]:
        """CI propagated from k_homing and k_death uncertainty."""
        if n0 <= 0:
            return [0.0, 0.0]
        ci_kh = self._p_ci["k_homing"]
        ci_kd = self._p_ci["k_death"]
        sigma_kh = (ci_kh[1] - ci_kh[0]) / (2 * 1.96)
        sigma_kd = (ci_kd[1] - ci_kd[0]) / (2 * 1.96)
        # Sensitivity: ∂delivered/∂k_homing, ∂delivered/∂k_death (numerical)
        k_total = self._p["k_homing"] + self._p["k_death"]
        # Analytical: delivered ≈ (k_homing / k_total) * n0 * stop_fraction
        frac = self._p["k_homing"] / k_total
        d_frac_dk_h = self._p["k_death"] / k_total ** 2
        d_frac_dk_d = -self._p["k_homing"] / k_total ** 2
        sens_h = n0 * self._p["stop_fraction"] * d_frac_dk_h
        sens_d = n0 * self._p["stop_fraction"] * d_frac_dk_d
        sigma = math.sqrt((sens_h * sigma_kh) ** 2 + (sens_d * sigma_kd) ** 2)
        return [round(delivered - 1.96 * sigma, 4), round(delivered + 1.96 * sigma, 4)]

    def _lag_ci95(self, lag_s: float) -> list[float]:
        """CI on lag from k_homing and k_death uncertainty."""
        if lag_s <= 0:
            return [0.0, 0.0]
        ci_kh = self._p_ci["k_homing"]
        ci_kd = self._p_ci["k_death"]
        sigma_kh = (ci_kh[1] - ci_kh[0]) / (2 * 1.96)
        sigma_kd = (ci_kd[1] - ci_kd[0]) / (2 * 1.96)
        k_total = self._p["k_homing"] + self._p["k_death"]
        # lag = -ln(1-stop_frac) / k_total
        d_lag_dk = math.log(1 - self._p["stop_fraction"]) / k_total ** 2  # negative
        sigma_lag = abs(d_lag_dk) * math.sqrt(sigma_kh ** 2 + sigma_kd ** 2) * 86_400
        return [round(lag_s - 1.96 * sigma_lag, 2), round(lag_s + 1.96 * sigma_lag, 2)]


# ------------------------------------------------------------------
# Entry point

def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(name)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="OISA — Blood Transit transfer model")
    parser.add_argument("--port", default="tcp://*:5011")
    parser.add_argument("--output-dir", default="logs/run/issl")
    args = parser.parse_args()

    model = BloodTransit(port=args.port, output_dir=Path(args.output_dir))
    model.run()


if __name__ == "__main__":
    main()
