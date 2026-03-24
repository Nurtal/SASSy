"""
OISSL / ISSL log loader — public analysis API.

Usage (matches examples in CLAUDE.md and README):

    from analysis.utils.oissl_parser import load_issl_series, load_oissl_series

    # Per-model ISSL logs
    records = load_issl_series("logs/COMP3/issl/bm_haematopoiesis/")
    hsc_counts = [r["continuous_state"][0]["count"] for r in records]

    # Orchestrator OISSL logs
    oissl = load_oissl_series("logs/COMP3/oissl/")
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator


def load_issl_series(directory: str | Path) -> list[dict]:
    """Load all ISSL checkpoint files from *directory*, sorted by sim_time_s.

    Files must follow the naming convention ``checkpoint_<sim_time_s>.json``.
    Returns a list of dicts (one per checkpoint), ascending by sim_time_s.
    """
    p = Path(directory)
    files = sorted(
        p.glob("checkpoint_*.json"),
        key=lambda f: int(f.stem.split("_")[1]),
    )
    return [json.loads(f.read_text()) for f in files]


def load_oissl_series(directory: str | Path) -> list[dict]:
    """Load all OISSL checkpoint files from *directory*, sorted by sim_time_s.

    Identical naming convention to ISSL — both functions are interchangeable
    for file loading purposes.
    """
    return load_issl_series(directory)


def iter_entity(
    records: list[dict],
    entity_id: str,
    section: str = "continuous_state",
) -> Iterator[tuple[float, dict]]:
    """Yield (sim_time_s, entity_dict) for every record containing *entity_id*.

    Useful for time-series extraction without writing list comprehensions:

        times, counts = zip(*[
            (t, e["count"]) for t, e in iter_entity(records, "CL:0000037")
        ])
    """
    for record in records:
        sim_time_s = record["envelope"]["sim_time_s"]
        for entity in record.get(section, []):
            if entity.get("entity_id") == entity_id or entity.get("signal_id") == entity_id:
                yield sim_time_s, entity
                break


def extract_export_flux(records: list[dict], signal_id: str) -> tuple[list[float], list[float]]:
    """Return (times_s, fluxes) for a named export signal across all records."""
    times, fluxes = [], []
    for record in records:
        t = record["envelope"]["sim_time_s"]
        for sig in record.get("export_signals", []):
            if sig["signal_id"] == signal_id:
                times.append(t)
                fluxes.append(sig["flux"])
                break
    return times, fluxes


def extract_composition_events(
    oissl_records: list[dict],
    event_type: str | None = None,
) -> list[dict]:
    """Flatten all composition_events from OISSL records, optionally filtered by type."""
    events = []
    for record in oissl_records:
        for ev in record.get("composition_events", []):
            if event_type is None or ev["type"] == event_type:
                events.append(ev)
    return events


def summary(records: list[dict]) -> dict:
    """Return a brief summary dict for a series of ISSL or OISSL records."""
    if not records:
        return {"n_checkpoints": 0}
    return {
        "n_checkpoints": len(records),
        "t_start_s":     records[0]["envelope"]["sim_time_s"],
        "t_end_s":       records[-1]["envelope"]["sim_time_s"],
        "model_id":      records[0]["envelope"].get("model_id",
                         records[0]["envelope"].get("run_id", "unknown")),
    }
