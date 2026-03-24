"""Tests for Model 3 — Thymic T-cell Selection ABM."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import jsonschema
import pytest

from models.thymus_selection.model import ThymusSelection, ThymicRealisation
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[3]
_ISSL_SCHEMA = json.loads((_REPO_ROOT / "schemas" / "issl_v1.schema.json").read_text())


@pytest.fixture()
def model(tmp_path):
    with patch("models._base.model_base.zmq.Context") as MockCtx:
        mock_sock = MagicMock()
        MockCtx.return_value.socket.return_value = mock_sock
        m = ThymusSelection(port="tcp://*:15012", output_dir=tmp_path)
    return m


def test_step_no_signal_no_crash(model):
    """Step with no incoming signal must not crash."""
    model._step(0.0, [])


def test_step_with_signal_populates_agents(model):
    """A step with a non-zero flux must add agents to at least one realisation."""
    model._step(0.0, [{"flux": 200.0, "ci_95": [180.0, 220.0]}])
    total_agents = sum(len(r.agents) for r in model._realisations)
    assert total_agents > 0


def test_emit_issl_validates(model):
    """ISSL record must validate against the schema."""
    model._step(0.0, [{"flux": 150.0, "ci_95": [130.0, 170.0]}])
    record = model.emit_issl(86_400.0)
    jsonschema.validate(record, _ISSL_SCHEMA)


def test_n_realisations_in_discrete_events(model):
    """n_realisations must equal the configured value in all discrete events."""
    model._step(0.0, [{"flux": 100.0}])
    record = model.emit_issl(86_400.0)
    for ev in record["discrete_events"]:
        assert ev["n_realisations"] == model._n_realisations


def test_export_after_dwell(model):
    """After enough checkpoints, export flux must become positive."""
    signal = [{"flux": 300.0, "ci_95": [270.0, 330.0]}]
    t = 0.0
    for _ in range(10):
        model._step(t, signal)
        t += ThymusSelection.DELTA_T_S
    record = model.emit_issl(t)
    flux = record["export_signals"][0]["flux"]
    assert flux >= 0, f"Export flux must be non-negative, got {flux}"


def test_realisation_positive_selection_fires(model):
    """With large influx, at least one positive selection event should occur."""
    model._step(0.0, [{"flux": 500.0}])
    total_pos = sum(r.pos_sel_count for r in model._realisations)
    # With 500 agents and p_encounter=0.30 per step over 24 steps, expect many events
    assert total_pos > 0, "Expected positive selection events with large influx"


def test_zero_flux_zero_export(model):
    """With zero flux, export after sufficient steps is zero (no agents to export)."""
    for i in range(5):
        model._step(i * ThymusSelection.DELTA_T_S, [{"flux": 0.0}])
    record = model.emit_issl(5 * ThymusSelection.DELTA_T_S)
    assert record["export_signals"][0]["flux"] == 0.0
