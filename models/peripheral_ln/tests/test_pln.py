"""Tests for Model 4 — Peripheral Lymph Node Dynamics."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import jsonschema
import pytest

from models.peripheral_ln.model import PeripheralLN

_REPO_ROOT = Path(__file__).resolve().parents[3]
_ISSL_SCHEMA = json.loads((_REPO_ROOT / "schemas" / "issl_v1.schema.json").read_text())


@pytest.fixture()
def model(tmp_path):
    with patch("models._base.model_base.zmq.Context") as MockCtx:
        mock_sock = MagicMock()
        MockCtx.return_value.socket.return_value = mock_sock
        m = PeripheralLN(port="tcp://*:15013", output_dir=tmp_path)
    return m


def test_step_no_crash(model):
    model._step(0.0, [])
    assert all(v >= 0 for v in model._state)


def test_emit_issl_validates(model):
    model._step(0.0, [])
    record = model.emit_issl(43_200.0)
    jsonschema.validate(record, _ISSL_SCHEMA)


def test_homeostasis_recovers_cd4(model):
    """Starting below set point, CD4 pool should grow toward S4."""
    model._state[0] = 50_000  # far below S4=200000
    initial_cd4 = float(model._state[0])
    for i in range(20):
        model._step(i * PeripheralLN.DELTA_T_S, [])
    assert model._state[0] > initial_cd4, "CD4 pool should grow toward homeostatic set point"


def test_import_signal_raises_pool(model):
    """Incoming thymic export signal should increase the pool over time."""
    # Start below set point so homeostatic proliferation and import both push up
    model._state[0] = 50_000
    model._state[1] = 30_000
    signal = [{"flux": 500.0}]
    for i in range(5):
        model._step(i * PeripheralLN.DELTA_T_S, signal)
    record = model.emit_issl(5 * PeripheralLN.DELTA_T_S)
    assert record["continuous_state"][0]["count"] > 50_000


def test_export_signals_empty_in_baseline(model):
    """export_signals must be empty for COMP-3 baseline configuration."""
    model._step(0.0, [])
    record = model.emit_issl(43_200.0)
    assert record["export_signals"] == []


def test_cd4_cd8_ratio_reasonable(model):
    """CD4/CD8 ratio at steady state should remain in biological range [1, 5]."""
    for i in range(40):
        model._step(i * PeripheralLN.DELTA_T_S, [])
    CD4, CD8 = model._state
    ratio = CD4 / CD8 if CD8 > 0 else float("inf")
    assert 1.0 <= ratio <= 5.0, f"CD4/CD8 ratio {ratio:.2f} outside expected range"
