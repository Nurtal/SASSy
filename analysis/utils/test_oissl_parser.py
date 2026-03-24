"""Tests for analysis.utils.oissl_parser."""

import json
import tempfile
from pathlib import Path

import pytest

from analysis.utils.oissl_parser import (
    extract_composition_events,
    extract_export_flux,
    iter_entity,
    load_issl_series,
    load_oissl_series,
    summary,
)


def _write_issl(directory: Path, sim_time_s: float, count: float, flux: float) -> None:
    record = {
        "envelope": {
            "model_id": "bm_haematopoiesis", "model_version": "3",
            "sim_time_s": sim_time_s, "schema_uri": "schemas/issl_v1.schema.json",
            "formalism": "ODE",
        },
        "continuous_state": [
            {"entity_id": "CL:0000037", "label": "HSC", "count": count,
             "unit": "cells", "ci_95": [count * 0.9, count * 1.1]},
        ],
        "discrete_events": [],
        "export_signals": [
            {"signal_id": "bm_haematopoiesis.progenitor_export",
             "entity_id": "CL:0002420", "flux": flux, "unit": "cells·day^-1",
             "lag_s": None, "ci_95": [flux * 0.9, flux * 1.1]},
        ],
        "internal_parameters": [],
        "watchdog": {"health_status": "ok", "ood_flag": False,
                     "divergence_score": None, "next_checkpoint_s": sim_time_s + 21600},
    }
    (directory / f"checkpoint_{int(sim_time_s)}.json").write_text(json.dumps(record))


@pytest.fixture()
def issl_dir(tmp_path):
    for t, count, flux in [(21600, 10000, 80), (43200, 10050, 82), (64800, 10100, 84)]:
        _write_issl(tmp_path, t, count, flux)
    return tmp_path


def test_load_issl_series_sorted(issl_dir):
    records = load_issl_series(issl_dir)
    assert len(records) == 3
    times = [r["envelope"]["sim_time_s"] for r in records]
    assert times == sorted(times)


def test_load_oissl_series_same_as_issl(issl_dir):
    assert load_oissl_series(issl_dir) == load_issl_series(issl_dir)


def test_iter_entity_yields_correct_entity(issl_dir):
    records = load_issl_series(issl_dir)
    results = list(iter_entity(records, "CL:0000037"))
    assert len(results) == 3
    times, entities = zip(*results)
    assert all(e["entity_id"] == "CL:0000037" for e in entities)


def test_extract_export_flux(issl_dir):
    records = load_issl_series(issl_dir)
    times, fluxes = extract_export_flux(records, "bm_haematopoiesis.progenitor_export")
    assert len(times) == 3
    assert fluxes == [80, 82, 84]


def test_summary(issl_dir):
    records = load_issl_series(issl_dir)
    s = summary(records)
    assert s["n_checkpoints"] == 3
    assert s["t_start_s"] == 21600
    assert s["t_end_s"] == 64800


def test_extract_composition_events_empty(issl_dir):
    # ISSL records have no composition_events (that's an OISSL field)
    records = load_issl_series(issl_dir)
    events = extract_composition_events(records)
    assert events == []
