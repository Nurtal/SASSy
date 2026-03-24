"""
Tests for Model 1 — Bone Marrow Haematopoiesis.

These tests instantiate the model *without* a live ZMQ socket by bypassing
the base class constructor (we test _step and emit_issl directly).
"""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import jsonschema
import pytest

from models.bm_haematopoiesis.model import BMHaematopoiesis

_REPO_ROOT = Path(__file__).resolve().parents[3]
_ISSL_SCHEMA = json.loads((_REPO_ROOT / "schemas" / "issl_v1.schema.json").read_text())


@pytest.fixture()
def model(tmp_path):
    """Return a BMHaematopoiesis instance with a mocked ZMQ socket."""
    with patch("models._base.model_base.zmq.Context") as MockCtx:
        mock_sock = MagicMock()
        MockCtx.return_value.socket.return_value = mock_sock
        m = BMHaematopoiesis(port="tcp://*:15010", output_dir=tmp_path)
    return m


def test_step_no_crash(model):
    """A single step must complete without error and keep state non-negative."""
    model._step(0.0, [])
    assert all(v >= 0 for v in model._state), "State variables must be non-negative."


def test_emit_issl_validates_against_schema(model):
    """emit_issl must produce a record that validates against issl_v1.schema.json."""
    model._step(0.0, [])
    record = model.emit_issl(21_600.0)
    # Raises jsonschema.ValidationError on failure
    jsonschema.validate(record, _ISSL_SCHEMA)


def test_export_flux_positive_after_30_days(model):
    """After 30 simulated days the DN1 export flux must be strictly positive."""
    seconds_per_day = 86_400
    checkpoint = BMHaematopoiesis.DELTA_T_S
    t = 0.0
    for _ in range(int(30 * seconds_per_day / checkpoint)):
        model._step(t, [])
        t += checkpoint
    record = model.emit_issl(t)
    flux = record["export_signals"][0]["flux"]
    assert flux > 0, f"Export flux should be positive after 30 days, got {flux}"


def test_watchdog_healthy_at_steady_state(model):
    """At steady state the watchdog should report health_status == 'ok'."""
    seconds_per_day = 86_400
    checkpoint = BMHaematopoiesis.DELTA_T_S
    t = 0.0
    for _ in range(int(60 * seconds_per_day / checkpoint)):
        model._step(t, [])
        t += checkpoint
    record = model.emit_issl(t)
    assert record["watchdog"]["health_status"] == "ok"


def test_ci95_bounds_contain_mean(model):
    """CI-95 lower bound must be ≤ count and upper bound must be ≥ count."""
    model._step(0.0, [])
    record = model.emit_issl(21_600.0)
    for entity in record["continuous_state"]:
        lo, hi = entity["ci_95"]
        assert lo <= entity["count"] <= hi, (
            f"{entity['label']}: count {entity['count']} not in CI [{lo}, {hi}]"
        )
