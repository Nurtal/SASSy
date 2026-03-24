"""
Tests for Component 2 — Scheduler.

We test the scheduler in isolation using mock ZMQ sockets that return
pre-baked ISSL records.  This verifies event ordering, signal delivery
timing, and OISSL emission without spawning real model processes.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

from orchestrator.components.causal_resolver import CausalResolver
from orchestrator.components.constraint_engine import ConstraintEngine
from orchestrator.components.ingestion import ISSLIngestion
from orchestrator.components.output_aggregator import OutputAggregator
from orchestrator.components.scheduler import Scheduler
from orchestrator.components.state_registry import StateRegistry
from orchestrator.components.transfer_dispatcher import TransferDispatcher
from orchestrator.components.watchdog import Watchdog

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _make_issl(model_id: str, sim_time_s: float, flux: float = 50.0) -> dict:
    return {
        "envelope": {
            "model_id": model_id, "model_version": "1",
            "sim_time_s": sim_time_s, "schema_uri": "schemas/issl_v1.schema.json",
            "formalism": "ODE",
        },
        "continuous_state": [
            {"entity_class": "obo", "entity_id": "CL:0000037",
             "label": "HSC", "count": 10000.0, "unit": "cells",
             "fitness": None, "surface_markers": [], "ci_95": [9000.0, 11000.0]},
        ],
        "discrete_events": [],
        "export_signals": [
            {"signal_id": f"{model_id}.progenitor_export",
             "entity_id": "CL:0002420", "flux": flux,
             "unit": "cells·day^-1", "lag_s": None, "ci_95": [40.0, 60.0]},
        ],
        "internal_parameters": [],
        "watchdog": {
            "health_status": "ok", "ood_flag": False,
            "divergence_score": None, "next_checkpoint_s": sim_time_s + 21600,
        },
    }


@pytest.fixture()
def components(tmp_path):
    """Build all 9 components with no live sockets or processes."""
    ingestion   = ISSLIngestion(output_dir=tmp_path)
    resolver    = CausalResolver(edges=[], transfer_model_ids=[])
    constraints = ConstraintEngine()
    registry    = StateRegistry()
    dispatcher  = TransferDispatcher(model_sockets={})
    aggregator  = OutputAggregator(
        output_dir=tmp_path, run_id="test",
        config_uri="configs/test.yaml", renderer_cfg={},
    )
    watchdog    = Watchdog()
    return dict(
        ingestion=ingestion, causal_resolver=resolver, constraint_engine=constraints,
        state_registry=registry, transfer_dispatcher=dispatcher,
        output_aggregator=aggregator, watchdog=watchdog,
    )


def _make_scheduler(components, tmp_path, global_clock, model_cfgs, sockets, edges=None):
    resolver = CausalResolver(edges=edges or [], transfer_model_ids=[])
    components["causal_resolver"] = resolver
    return Scheduler(
        global_clock=global_clock,
        model_cfgs=model_cfgs,
        model_sockets=sockets,
        **components,
    )


def test_scheduler_runs_single_model_two_checkpoints(components, tmp_path):
    """Scheduler must fire model checkpoints and emit OISSL at global intervals."""
    model_cfg = [{"id": "bm", "delta_t_s": 43200.0, "issl_port": "tcp://localhost:9999"}]
    global_clock = {"start_s": 0, "end_s": 86400, "checkpoint_interval_s": 86400}

    # Mock socket: alternates ack / issl_record on recv_json
    sock = MagicMock()
    responses = []
    # 2 steps (t=43200, t=86400) each needs: ack for step, issl for emit
    for t in [43200.0, 86400.0]:
        responses.append({"status": "ok"})   # step ack
        responses.append(_make_issl("bm", t))  # emit record
    sock.recv_json.side_effect = responses

    sched = _make_scheduler(components, tmp_path, global_clock, model_cfg,
                             sockets={"bm": sock})
    sched.run()

    # Should have produced 1 OISSL file at t=86400
    oissl_files = list((tmp_path / "oissl").glob("checkpoint_*.json"))
    assert len(oissl_files) == 1
    assert "86400" in oissl_files[0].name


def test_scheduler_emits_oissl_at_global_interval(components, tmp_path):
    """With checkpoint_interval_s=43200, two OISSL files should exist after 86400 s."""
    model_cfg = [{"id": "bm", "delta_t_s": 43200.0, "issl_port": "tcp://localhost:9999"}]
    global_clock = {"start_s": 0, "end_s": 86400, "checkpoint_interval_s": 43200}

    sock = MagicMock()
    responses = []
    for t in [43200.0, 86400.0]:
        responses.append({"status": "ok"})
        responses.append(_make_issl("bm", t))
    sock.recv_json.side_effect = responses

    sched = _make_scheduler(components, tmp_path, global_clock, model_cfg,
                             sockets={"bm": sock})
    sched.run()

    oissl_files = list((tmp_path / "oissl").glob("checkpoint_*.json"))
    assert len(oissl_files) == 2


def test_state_registry_updated_after_model_step(components, tmp_path):
    """StateRegistry must contain the model's entities after one checkpoint."""
    model_cfg = [{"id": "bm", "delta_t_s": 86400.0, "issl_port": "tcp://localhost:9999"}]
    global_clock = {"start_s": 0, "end_s": 86400, "checkpoint_interval_s": 86400}

    sock = MagicMock()
    sock.recv_json.side_effect = [
        {"status": "ok"},
        _make_issl("bm", 86400.0),
    ]

    sched = _make_scheduler(components, tmp_path, global_clock, model_cfg,
                             sockets={"bm": sock})
    sched.run()

    snap = components["state_registry"].get("CL:0000037")
    assert snap is not None
    assert snap.count == 10000.0
