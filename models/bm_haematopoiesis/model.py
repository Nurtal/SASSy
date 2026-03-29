#!/usr/bin/env python3
"""
Model 1 — Bone Marrow Haematopoiesis (ODE, v7 — 5-compartment)

State vector: [HSC, MPP, LMPP, CLP, DN1]

ODEs (all rates in day⁻¹; converted from seconds before integration):

  dHSC/dt  = r_self × HSC × (1 − HSC/K_niche) − (d_HSC_MPP + d_apop) × HSC
  dMPP/dt  = d_HSC_MPP × HSC + r_MPP × MPP − (d_MPP_LMPP + d_MPP_death) × MPP
  dLMPP/dt = f_lymphoid × d_MPP_LMPP × MPP − (d_LMPP_CLP + d_LMPP_death) × LMPP
  dCLP/dt  = d_LMPP_CLP × LMPP − (alpha_T + d_CLP_other + d_CLP_death) × CLP
  dDN1/dt  = alpha_T × CLP − export_rate × DN1

Export flux (cells/day): export_rate × DN1  ≈ 59.55 at equilibrium
  (within murine literature range 10–100 cells·day⁻¹; Bhandoola et al. 2007)

Analytical steady-state (computed from parameters; v7 exact floats):
  HSC  =  9075.0000   (consistent with Bhatt 2016: ~10,000 LSK CD150+CD48-)
  MPP  = 27225.0000   (consistent with Wilson 2008: ~15,000–50,000 LSK CD150-)
  LMPP = 11910.9375   (consistent with Adolfsson 2005: ~10,000–20,000 LSK Flt3+)
  CLP  = 29777.3438   (consistent with Rodrigues 2005: ~20,000–30,000 Lin-IL7Rα+Flt3+)
  DN1  =  1191.0938   (ETP pool in BM before blood egress)
  flux =    59.5547   cells·day⁻¹

v6 → v7 change: CI-95 replaced by Monte Carlo parameter propagation.
  Why: the near-balanced MPP regime (λ_MPP = 0.001 day⁻¹) makes first-order
  linear CI propagation invalid.  The true sensitivity of MPP to r_MPP is
  S = r_MPP/λ_MPP ≈ 199, which the former S=1 formula underestimated by 2 orders
  of magnitude.  Moreover ~48% of the parameter 95%-CI space gives λ_MPP ≤ 0
  (unstable/divergent MPP), so no finite analytical CI exists for MPP or any
  downstream compartment.  _compute_mc_ci95() samples 2 000 parameter sets
  from N(μ, σ) (σ = CI-width/3.92), retains stable draws (λ_MPP > 0), and
  reports the 2.5th/97.5th percentiles.  Stability fraction (~51%) is stored in
  the ISSL watchdog.

  MC CI at nominal parameters (2 000 stable draws):
    HSC  : [ 7 235,  10 894]
    MPP  : [   570,  43 603]   (heavily right-skewed; median ~2 000)
    LMPP : [   269,  19 350]
    CLP  : [   650,  48 557]
    DN1  : [    22,   2 006]
    flux : [  1.1,   101.3]  cells·day⁻¹  (spans full literature range)

Transit amplification: r_MPP ≈ d_MPP_LMPP + d_MPP_death (near-balanced).
  λ_MPP = d_MPP_LMPP + d_MPP_death − r_MPP = 0.001 day⁻¹  (stability: λ_MPP > 0 required)

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
    MODEL_VERSION = "7"
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

        # Initialise at analytical steady-state (exact floats, not rounded integers).
        self._state = self._compute_steady_state()
        self._sim_time_s: float = 0.0

        # Pre-compute Monte Carlo CI-95 once (parameter uncertainty propagation).
        # Stored as a dict: {"HSC":[lo,hi], "MPP":[lo,hi], ..., "flux":[lo,hi],
        #                    "n_stable": int, "n_total": int}
        self._ci95 = self._compute_mc_ci95()

    # ------------------------------------------------------------------
    # Analytical steady-state initialisation

    def _compute_steady_state(self) -> np.ndarray:
        """Return [HSC, MPP, LMPP, CLP, DN1] at analytical fixed point (dX/dt = 0).

        Closed-form cascade solution:
          HSC*  = K_niche × (1 − (d_HSC_MPP + d_apop) / r_self)
          MPP*  = d_HSC_MPP × HSC* / λ_MPP           where λ_MPP = d_MPP_LMPP + d_MPP_death − r_MPP
          LMPP* = f_lymphoid × d_MPP_LMPP × MPP* / (d_LMPP_CLP + d_LMPP_death)
          CLP*  = d_LMPP_CLP × LMPP* / (alpha_T + d_CLP_other + d_CLP_death)
          DN1*  = alpha_T × CLP* / export_rate

        Raises ValueError if λ_MPP ≤ 0 (unstable / unbounded MPP pool).
        """
        p = self._p
        HSC  = p["K_niche"] * (1.0 - (p["d_HSC_MPP"] + p["d_apop"]) / p["r_self"])
        lam  = p["d_MPP_LMPP"] + p["d_MPP_death"] - p["r_MPP"]
        if lam <= 0.0:
            raise ValueError(
                f"MPP net loss λ_MPP = {lam:.6g} day⁻¹ ≤ 0 — "
                "system is unstable; check r_MPP < d_MPP_LMPP + d_MPP_death"
            )
        MPP  = p["d_HSC_MPP"] * HSC / lam
        LMPP = p["f_lymphoid"] * p["d_MPP_LMPP"] * MPP / (p["d_LMPP_CLP"] + p["d_LMPP_death"])
        CLP  = p["d_LMPP_CLP"] * LMPP / (p["alpha_T"] + p["d_CLP_other"] + p["d_CLP_death"])
        DN1  = p["alpha_T"] * CLP / p["export_rate"]
        return np.array([HSC, MPP, LMPP, CLP, DN1])

    def _compute_mc_ci95(self, n_samples: int = 2_000, seed: int = 42) -> dict:
        """Monte Carlo CI-95 via vectorised parameter sampling.

        Samples all parameters jointly from N(μ, σ) where σ = (ci_hi−ci_lo)/3.92.
        Draws are rejected when λ_MPP = d_MPP_LMPP + d_MPP_death − r_MPP ≤ 0
        (unstable regime) or when any parameter is non-positive.  The 2.5th and
        97.5th percentiles of the retained draws give the CI.

        Returns
        -------
        dict with keys "HSC","MPP","LMPP","CLP","DN1","flux" → [lo, hi]
             and "n_stable" (int), "n_total" (int).
        """
        rng  = np.random.default_rng(seed)
        keys = list(self._p.keys())
        mu   = np.array([self._p[k]        for k in keys])
        cis  = np.array([self._p_ci[k]     for k in keys])
        sig  = (cis[:, 1] - cis[:, 0]) / (2.0 * 1.96)

        names  = ("HSC", "MPP", "LMPP", "CLP", "DN1", "flux")
        store  = {k: [] for k in names}
        n_total = 0
        chunk   = 500

        while len(store["HSC"]) < n_samples and n_total < n_samples * 60:
            batch = rng.normal(mu, sig, size=(chunk, len(keys)))
            n_total += chunk
            pp = {k: batch[:, i] for i, k in enumerate(keys)}

            lam   = pp["d_MPP_LMPP"] + pp["d_MPP_death"] - pp["r_MPP"]
            HSC_v = pp["K_niche"] * (1.0 - (pp["d_HSC_MPP"] + pp["d_apop"]) / pp["r_self"])
            valid = (lam > 0) & np.all(batch > 0, axis=1) & (HSC_v > 0)

            lam   = lam[valid]
            HSC_v = HSC_v[valid]
            pv    = {k: pp[k][valid] for k in keys}

            MPP_v  = pv["d_HSC_MPP"] * HSC_v / lam
            LMPP_v = pv["f_lymphoid"] * pv["d_MPP_LMPP"] * MPP_v / (pv["d_LMPP_CLP"] + pv["d_LMPP_death"])
            CLP_v  = pv["d_LMPP_CLP"] * LMPP_v / (pv["alpha_T"] + pv["d_CLP_other"] + pv["d_CLP_death"])
            DN1_v  = pv["alpha_T"] * CLP_v / pv["export_rate"]
            flux_v = pv["export_rate"] * DN1_v

            for name, arr in zip(names, [HSC_v, MPP_v, LMPP_v, CLP_v, DN1_v, flux_v]):
                store[name].extend(arr.tolist())

        result: dict = {
            "n_stable": len(store["HSC"]),
            "n_total":  n_total,
        }
        for name in names:
            arr = np.array(store[name][:n_samples])
            result[name] = [round(float(np.percentile(arr, 2.5)),  2),
                            round(float(np.percentile(arr, 97.5)), 2)]
        return result

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

        # CI-95 — Monte Carlo parameter propagation (pre-computed in __init__).
        # See _compute_mc_ci95 and module docstring for rationale.
        ci_hsc    = self._ci95["HSC"]
        ci_mpp    = self._ci95["MPP"]
        ci_lmpp   = self._ci95["LMPP"]
        ci_clp    = self._ci95["CLP"]
        ci_dn1    = self._ci95["DN1"]
        ci_export = self._ci95["flux"]

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
            "watchdog": {
                **self._make_watchdog(
                    sim_time_s,
                    ood_flag=self._is_ood(HSC, MPP, LMPP, CLP, DN1),
                    divergence_score=self._divergence(HSC),
                ),
                "ci_method": "monte_carlo",
                "ci_n_stable": self._ci95["n_stable"],
                "ci_n_total":  self._ci95["n_total"],
                "ci_stability_fraction": round(
                    self._ci95["n_stable"] / max(self._ci95["n_total"], 1), 4
                ),
                "lambda_mpp": round(
                    p["d_MPP_LMPP"] + p["d_MPP_death"] - p["r_MPP"], 6
                ),
            },
        }

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
        """Out-of-domain flag based on v5 equilibrium bounds."""
        return bool(
            HSC  < 500    or HSC  > 50_000
            or MPP  < 100    or MPP  > 500_000
            or LMPP < 100    or LMPP > 200_000
            or CLP  < 100    or CLP  > 300_000
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
    parser = argparse.ArgumentParser(description="OISA — BM Haematopoiesis model (v7)")
    parser.add_argument("--port",       default="tcp://*:5010",
                        help="ZMQ REP bind address")
    parser.add_argument("--output-dir", default="logs/run/issl",
                        help="Base directory for ISSL checkpoint files")
    args = parser.parse_args()

    model = BMHaematopoiesis(port=args.port, output_dir=Path(args.output_dir))
    model.run()


if __name__ == "__main__":
    main()
