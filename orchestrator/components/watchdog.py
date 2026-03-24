"""
Component 9 — Watchdog Monitor

Inspects the watchdog section of each emitted ISSL record.
Returns WatchdogAlert objects that are added to OISSL composition_events.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class WatchdogAlert:
    model_id: str
    health_status: str
    ood_flag: bool
    divergence_score: float | None
    sim_time_s: float
    message: str


class Watchdog:
    """Health monitor for model processes."""

    def __init__(self) -> None:
        self._alerts: list[WatchdogAlert] = []

    def check(self, model_id: str, record: dict, sim_time_s: float) -> WatchdogAlert | None:
        """Inspect *record*'s watchdog section. Returns an alert or None."""
        w = record.get("watchdog", {})
        health = w.get("health_status", "ok")
        ood = bool(w.get("ood_flag", False))
        divergence = w.get("divergence_score")

        alert: WatchdogAlert | None = None

        if health != "ok":
            msg = f"Model {model_id} health_status={health}"
            if divergence is not None:
                msg += f" (divergence={divergence:.4f})"
            alert = WatchdogAlert(
                model_id=model_id,
                health_status=health,
                ood_flag=ood,
                divergence_score=divergence,
                sim_time_s=sim_time_s,
                message=msg,
            )
            logger.warning("WATCHDOG: %s", msg)

        elif ood:
            msg = f"Model {model_id} out-of-distribution (divergence={divergence})"
            alert = WatchdogAlert(
                model_id=model_id,
                health_status="ok",
                ood_flag=True,
                divergence_score=divergence,
                sim_time_s=sim_time_s,
                message=msg,
            )
            logger.warning("WATCHDOG OOD: %s", msg)

        if alert:
            self._alerts.append(alert)

        return alert

    def alerts_since(self, sim_time_s: float) -> list[WatchdogAlert]:
        return [a for a in self._alerts if a.sim_time_s >= sim_time_s]

    def to_composition_events(self, alerts: list[WatchdogAlert]) -> list[dict]:
        return [
            {
                "type":         "watchdog_alert",
                "sim_time_s":   a.sim_time_s,
                "source_model": a.model_id,
                "message":      a.message,
                "value":        a.divergence_score,
                "severity":     "error" if a.health_status == "failed" else "warning",
            }
            for a in alerts
        ]

    def report(self) -> None:
        """Log a summary of all alerts at end of run."""
        if not self._alerts:
            logger.info("Watchdog: no alerts during run.")
            return
        logger.warning("Watchdog: %d alert(s) during run:", len(self._alerts))
        for a in self._alerts:
            logger.warning("  t=%-10.0f  %-30s  %s", a.sim_time_s, a.model_id, a.message)
