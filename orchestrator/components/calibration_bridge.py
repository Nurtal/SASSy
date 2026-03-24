"""
Component 7 — Calibration Bridge

Provides externally-calibrated signal values (e.g. from EHR time-series,
experimental data, or basal constants) to models that request them.

In the COMP-3 baseline configuration all signals return None (basal values
are baked into model parameters).  Future implementations can inject
patient-specific EHR time-series by subclassing or configuring this bridge.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class CalibrationBridge:
    """External signal provider.  Default: all queries return None (use model defaults)."""

    def __init__(self, data_path: str | Path | None = None) -> None:
        self._series: dict[str, pd.Series] = {}

        if data_path is not None:
            self._load(Path(data_path))

    def _load(self, path: Path) -> None:
        """Load a CSV with columns: signal_name, sim_time_s, value."""
        if not path.exists():
            logger.warning("CalibrationBridge: data file not found at %s — using defaults", path)
            return
        df = pd.read_csv(path)
        for signal_name, grp in df.groupby("signal_name"):
            self._series[signal_name] = (
                grp.set_index("sim_time_s")["value"].sort_index()
            )
        logger.info("CalibrationBridge: loaded %d signal(s) from %s",
                    len(self._series), path)

    def get_signal(self, signal_name: str, sim_time_s: float) -> float | None:
        """Return the calibrated value for *signal_name* at *sim_time_s*, or None.

        If a time-series is loaded, the nearest-index value is returned.
        Returns None when the signal is not configured (model uses its own default).
        """
        if signal_name not in self._series:
            return None
        series = self._series[signal_name]
        # Nearest-index lookup
        idx = series.index.get_indexer([sim_time_s], method="nearest")[0]
        return float(series.iloc[idx])

    def inject_as_signal(self, signal_name: str, sim_time_s: float) -> dict | None:
        """Return a signal dict for injection into a model's step(), or None."""
        value = self.get_signal(signal_name, sim_time_s)
        if value is None:
            return None
        return {
            "signal_id": signal_name,
            "entity_id": signal_name,
            "flux":      value,
            "unit":      "calibrated",
        }
