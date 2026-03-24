"""Tests for Model 2 — Blood Transit Kinetics (transfer model)."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import jsonschema
import pytest

from models.blood_transit.model import BloodTransit

_REPO_ROOT = Path(__file__).resolve().parents[3]
_ISSL_SCHEMA = json.loads((_REPO_ROOT / "schemas" / "issl_v1.schema.json").read_text())


@pytest.fixture()
def model(tmp_path):
    with patch("models._base.model_base.zmq.Context") as MockCtx:
        mock_sock = MagicMock()
        MockCtx.return_value.socket.return_value = mock_sock
        m = BloodTransit(port="tcp://*:15011", output_dir=tmp_path)
    return m


def test_step_with_signal_computes_delivery(model):
    """Transit model must produce positive delivered cells from a positive flux."""
    signal = {"signal_id": "bm_haematopoiesis.progenitor_export",
               "flux": 100.0, "unit": "cells·day^-1"}
    model._step(0.0, [signal])
    assert model._last_cells_delivered > 0
    assert model._last_lag_s > 0


def test_lag_is_in_seconds(model):
    """lag_s must be in seconds (reasonable transit: 12 h – 7 days)."""
    signal = {"flux": 200.0}
    model._step(0.0, [signal])
    assert 3_600 <= model._last_lag_s <= 7 * 86_400, (
        f"lag_s={model._last_lag_s} outside expected 1 h – 7 days range"
    )


def test_conservation_cells_delivered_leq_n0(model):
    """Cells delivered can never exceed the initial pool size."""
    n0 = 500.0
    model._step(0.0, [{"flux": n0}])
    assert model._last_cells_delivered <= n0 + 1e-6  # tolerance for numerics


def test_emit_issl_validates(model):
    """ISSL record from emit_issl must validate against the schema."""
    model._step(0.0, [{"flux": 80.0}])
    record = model.emit_issl(0.0)
    jsonschema.validate(record, _ISSL_SCHEMA)


def test_lag_s_present_in_export_signals(model):
    """export_signals[0].lag_s must be populated (non-null) after a step with flux."""
    model._step(0.0, [{"flux": 50.0}])
    record = model.emit_issl(0.0)
    lag = record["export_signals"][0]["lag_s"]
    assert lag is not None and lag > 0


def test_zero_flux_returns_zero_delivery(model):
    """Zero input flux must produce zero delivered cells and zero lag."""
    model._step(0.0, [{"flux": 0.0}])
    assert model._last_cells_delivered == 0.0
    assert model._last_lag_s == 0.0
