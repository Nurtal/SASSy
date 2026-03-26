#!/usr/bin/env python3
"""
Model 1 — Bone Marrow Haematopoiesis (ODE, v4 — 5-compartment)

State vector: [HSC, MPP, LMPP, CLP, DN1]

ODEs (all rates in day⁻¹; converted from seconds before integration):

  dHSC/dt  = r_self × HSC × (1 − HSC/K_niche) − (d_HSC_MPP + d_apop) × HSC
  dMPP/dt  = d_HSC_MPP × HSC + r_MPP × MPP − (d_MPP_LMPP + d_MPP_death) × MPP
  dLMPP/dt = f_lymphoid × d_MPP_LMPP × MPP − (d_LMPP_CLP + d_LMPP_death) × LMPP
  dCLP/dt  = d_LMPP_CLP × LMPP − (alpha_T + d_CLP_other + d_CLP_death) × CLP
  dDN1/dt  = alpha_T × CLP − export_rate × DN1

Export flux (cells/day): export_rate × DN1  ≈ 27.2 at equilibrium
  (4.3× the v3 value; within murine literature range 10–100 cells·day⁻¹)

Key equilibrium values (v4):
  HSC = 9075, MPP = 1361, LMPP = 595, CLP = 340, DN1 = 181

Transit amplification: r_MPP > 0 allows MPP export flux to exceed HSC input flux.
  λ_MPP = d_MPP_LMPP + d_MPP_death − r_MPP = 0.020 day⁻¹  (stability: λ_MPP > 0 required)

CI-95 note for MPP: the near-balanced regime means true Jacobian sensitivity
  coefficients S(r_MPP) = r_MPP/λ_MPP ≈ 10.  _mpp_ci95 uses the standard
  relative-σ formula (S = 1) as a conservative lower bound on uncertainty.

Orchestrator compatibility:
  export_signals[0].signal_id = "bm_haematopoiesis.progenitor_export"  (unchanged)
  export_signals[0].entity_id = "CL:0002420"                           (unchanged)

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

# OBO Cell Ontology entity IDs
_ID_HSC  = "CL:0000037"   # haematopoietic stem cell
_ID_MPP  = "CL:0000837"   # haematopoietic multipotent progenitor cell
_ID_LMPP = "CL:0000838"   # lymphoid-primed multipotent progenitor
_ID_CLP  = "CL:0000051"   # common lymphoid progenitor
_ID_DN1  = "CL:0002420"   # early T-lineage progenitor (DN1/ETP)


def _load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


class BMHaematopoiesis(ModelBase):
    """Bone marrow 5-compartment ODE model (HSC → MPP → LMPP → CLP → DN1)."""

    MODEL_ID      = "bm_haematopoiesis"
    MODEL_VERSION = "4"
    DELTA_T_S     = 21_600.0   # 6 hours in seconds

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
        self._p      = {k: v["value"] for k, v in cfg["parameters"].items()}
        self._p_ci   = {k: v["ci_95"]  for k, v in cfg["parameters"].items()}
        self._p_meta = cfg["parameters"]

        ic = cfg["initial_conditions"]
        # State: [HSC, MPP, LMPP, CLP, DN1]
        self._state = np.array(
            [ic["HSC"], ic["MPP"], ic["LMPP"], ic["CLP"], ic["DN1"]],
            dtype=float,
        )
        self._sim_time_s: float = 0.0

    # ------------------------------------------------------------------
    # ODE system

    def _ode(self, t: float, y: np.ndarray) -> np.ndarray:
        HSC, MPP, LMPP, CLP, DN1 = y
        p = self._p
        dHSC  = (p["r_self"] * HSC * (1.0 - HSC / p["K_niche"])
                 - (p["d_HSC_MPP"] + p["d_apop"]) * HSC)
        dMPP  = (p["d_HSC_MPP"] * HSC
                 + p["r_MPP"] * MPP
                 - (p["d_MPP_LMPP"] + p["d_MPP_death"]) * MPP)
        dLMPP = (p["f_lymphoid"] * p["d_MPP_LMPP"] * MPP
                 - (p["d_LMPP_CLP"] + p["d_LMPP_death"]) * LMPP)
        dCLP  = (p["d_LMPP_CLP"] * LMPP
                 - (p["alpha_T"] + p["d_CLP_other"] + p["d_CLP_death"]) * CLP)
        dDN1  = p["alpha_T"] * CLP - p["export_rate"] * DN1
        return np.array([dHSC, dMPP, dLMPP, dCLP, dDN1])

    # ------------------------------------------------------------------
    # Abstract interface

    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        """Advance the ODE by one DELTA_T_S starting from *sim_time_s*."""
        _ = signals   # niche cytokine signal — unused in baseline

        t_start_d = sim_time_s / 86_400
        t_end_d   = (sim_time_s + self.DELTA_T_S) / 86_400

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

        self._state    = np.maximum(sol.y[:, -1], 0.0)
        self._sim_time_s = sim_time_s + self.DELTA_T_S

    def emit_issl(self, sim_time_s: float) -> dict:
        HSC, MPP, LMPP, CLP, DN1 = self._state
        p = self._p

        export_flux = p["export_rate"] * DN1   # cells/day

        # CI-95 — first-order relative-σ propagation (see module docstring)
        ci_export = self._flux_ci95(DN1, export_flux)
        ci_hsc    = self._state_ci95(HSC,  "r_self", "K_niche")
        ci_mpp    = self._mpp_ci95(MPP)
        ci_lmpp   = self._state_ci95(LMPP, "f_lymphoid", "d_MPP_LMPP",
                                     "d_LMPP_CLP", "d_LMPP_death")
        ci_clp    = self._state_ci95(CLP,  "d_LMPP_CLP", "alpha_T",
                                     "d_CLP_other", "d_CLP_death")
        ci_dn1    = self._state_ci95(DN1,  "alpha_T", "export_rate")

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
                    "entity_class":    "obo",
                    "entity_id":       _ID_HSC,
                    "label":           "Haematopoietic stem cell",
                    "count":           round(HSC, 2),
                    "unit":            "cells",
                    "fitness":         None,
                    "surface_markers": ["CD34+", "CD117+", "CD150+", "CD48-", "Sca1+"],
                    "ci_95":           ci_hsc,
                },
                {
                    "entity_class":    "obo",
                    "entity_id":       _ID_MPP,
                    "label":           "Multipotent progenitor",
                    "count":           round(MPP, 2),
                    "unit":            "cells",
                    "fitness":         None,
                    "surface_markers": ["CD34+", "CD117+", "CD150-", "CD48+", "Flt3-"],
                    "ci_95":           ci_mpp,
                },
                {
                    "entity_class":    "obo",
                    "entity_id":       _ID_LMPP,
                    "label":           "Lymphoid-primed multipotent progenitor",
                    "count":           round(LMPP, 2),
                    "unit":            "cells",
                    "fitness":         None,
                    "surface_markers": ["CD34+", "CD117lo", "Flt3+", "Sca1+", "CD127-"],
                    "ci_95":           ci_lmpp,
                },
                {
                    "entity_class":    "obo",
                    "entity_id":       _ID_CLP,
                    "label":           "Common lymphoid progenitor",
                    "count":           round(CLP, 2),
                    "unit":            "cells",
                    "fitness":         None,
                    "surface_markers": ["CD34+", "CD127+", "CD117lo", "Flt3+"],
                    "ci_95":           ci_clp,
                },
                {
                    "entity_class":    "obo",
                    "entity_id":       _ID_DN1,
                    "label":           "Early T-lineage progenitor (DN1/ETP)",
                    "count":           round(DN1, 2),
                    "unit":            "cells",
                    "fitness":         None,
                    "surface_markers": ["CD34+", "CD117+", "CD44+", "CD25-"],
                    "ci_95":           ci_dn1,
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
                ood_flag=self._is_ood(HSC, MPP, LMPP, CLP, DN1),
                divergence_score=self._divergence(HSC),
            ),
        }

    # ------------------------------------------------------------------
    # CI-95 helpers

    def _state_ci95(self, value: float, *param_keys: str) -> list[float]:
        """Approximate CI-95 using first-order relative-σ (diagonal covariance).

        For each parameter p_i with nominal value p_i0 and CI [lo, hi]:
            σ_rel_i = (hi − lo) / (2 × 1.96 × p_i0)
        Combined relative σ:
            σ_rel = sqrt(Σ σ_rel_i²)
        Absolute σ on the state:
            σ_abs = σ_rel × value

        All sensitivity coefficients assumed = 1 (first-order linear approx).
        For MPP in the near-balanced regime, use _mpp_ci95 instead.
        """
        variance_rel = 0.0
        for k in param_keys:
            ci        = self._p_ci[k]
            param_val = self._p[k]
            sigma_rel = (ci[1] - ci[0]) / (2.0 * 1.96 * param_val)
            variance_rel += sigma_rel ** 2
        sigma_abs = math.sqrt(variance_rel) * value
        return [round(value - 1.96 * sigma_abs, 2), round(value + 1.96 * sigma_abs, 2)]

    def _mpp_ci95(self, mpp: float) -> list[float]:
        """CI-95 for MPP — conservative lower bound for near-balanced regime.

        The true Jacobian sensitivity coefficients are:
            S(r_MPP)      = r_MPP      / λ_MPP ≈ 0.20/0.02 = 10
            S(d_MPP_LMPP) = d_MPP_LMPP / λ_MPP ≈ 0.15/0.02 =  7.5
            S(d_MPP_death)= d_MPP_death / λ_MPP ≈ 0.07/0.02 =  3.5

        Using S = 10 gives CI that reaches negative values.  We use S = 1
        (standard _state_ci95) as a conservative lower bound.  The result is
        narrower than the true posterior, so treat MPP CI as indicative only.
        """
        return self._state_ci95(mpp, "r_MPP", "d_MPP_LMPP", "d_MPP_death")

    def _flux_ci95(self, dn1: float, flux: float) -> list[float]:
        """CI-95 for export flux = export_rate × DN1, from export_rate uncertainty."""
        p    = self._p
        ci   = self._p_ci["export_rate"]
        σ_er = (ci[1] - ci[0]) / (2.0 * 1.96)
        σ_f  = dn1 * σ_er
        return [round(flux - 1.96 * σ_f, 4), round(flux + 1.96 * σ_f, 4)]

    # ------------------------------------------------------------------
    # Watchdog helpers

    def _is_ood(
        self,
        HSC: float,
        MPP: float,
        LMPP: float,
        CLP: float,
        DN1: float,
    ) -> bool:
        """Out-of-domain flag: HSC niche bounds + non-negativity for all others."""
        return bool(
            HSC  < 500 or HSC > 50_000
            or MPP  < 0
            or LMPP < 0
            or CLP  < 0
            or DN1  < 0
        )

    def _divergence(self, HSC: float) -> float:
        """Scalar divergence score: relative deviation of HSC from niche set point."""
        return abs(HSC - self._p["K_niche"]) / self._p["K_niche"]


# ------------------------------------------------------------------
# Entry point

def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )
    parser = argparse.ArgumentParser(description="OISA — BM Haematopoiesis model (v4)")
    parser.add_argument("--port",       default="tcp://*:5010",
                        help="ZMQ REP bind address")
    parser.add_argument("--output-dir", default="logs/run/issl",
                        help="Base directory for ISSL checkpoint files")
    args = parser.parse_args()

    model = BMHaematopoiesis(port=args.port, output_dir=Path(args.output_dir))
    model.run()


if __name__ == "__main__":
    main()
