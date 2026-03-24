"""Tests for Component 3 — CausalResolver."""

import pytest
from orchestrator.components.causal_resolver import CausalResolver


_EDGES = [
    {
        "source": "bm_haematopoiesis",
        "signal_id": "bm_haematopoiesis.progenitor_export",
        "target": "thymus_selection",
        "lag": "constant:86400",
        "activation_threshold": 0.0,
    },
    {
        "source": "thymus_selection",
        "signal_id": "thymus_selection.naive_T_export",
        "target": "peripheral_ln",
        "lag": "constant:172800",
        "activation_threshold": 0.0,
    },
]

_ISSL_BM = {
    "export_signals": [
        {
            "signal_id": "bm_haematopoiesis.progenitor_export",
            "entity_id": "CL:0002420",
            "flux": 80.0,
            "unit": "cells·day^-1",
            "lag_s": None,
            "ci_95": [70.0, 90.0],
        }
    ]
}


def test_dag_builds_without_error():
    resolver = CausalResolver(_EDGES, transfer_model_ids=[])
    assert resolver is not None


def test_cycle_detection():
    cyclic_edges = [
        {"source": "A", "signal_id": "A.sig", "target": "B", "lag": "constant:0"},
        {"source": "B", "signal_id": "B.sig", "target": "A", "lag": "constant:0"},
    ]
    with pytest.raises(ValueError, match="cycle"):
        CausalResolver(cyclic_edges, transfer_model_ids=[])


def test_get_edges_from_source():
    resolver = CausalResolver(_EDGES, transfer_model_ids=[])
    edges = resolver.get_edges_from("bm_haematopoiesis")
    assert len(edges) == 1
    assert edges[0]["target"] == "thymus_selection"


def test_route_constant_lag():
    resolver = CausalResolver(_EDGES, transfer_model_ids=[])
    routed = resolver.route("bm_haematopoiesis", _ISSL_BM, sim_time_s=0.0,
                             transfer_issl_records={})
    assert len(routed) == 1
    rs = routed[0]
    assert rs.target_model == "thymus_selection"
    assert rs.lag_s == 86400.0
    assert rs.deliver_at_s == 86400.0


def test_route_below_threshold_not_routed():
    edges = [{
        "source": "bm_haematopoiesis",
        "signal_id": "bm_haematopoiesis.progenitor_export",
        "target": "thymus_selection",
        "lag": "constant:86400",
        "activation_threshold": 100.0,  # flux=80 < 100 → not routed
    }]
    resolver = CausalResolver(edges, transfer_model_ids=[])
    routed = resolver.route("bm_haematopoiesis", _ISSL_BM, 0.0, {})
    assert len(routed) == 0


def test_route_model_lag_uses_transfer_record():
    edges = [{
        "source": "bm_haematopoiesis",
        "signal_id": "bm_haematopoiesis.progenitor_export",
        "target": "thymus_selection",
        "lag": "model:blood_transit",
        "activation_threshold": 0.0,
    }]
    resolver = CausalResolver(edges, transfer_model_ids=["blood_transit"])
    transfer_records = {
        "blood_transit": {
            "export_signals": [{"lag_s": 57600.0, "flux": 75.0}]
        }
    }
    routed = resolver.route("bm_haematopoiesis", _ISSL_BM, 0.0, transfer_records)
    assert len(routed) == 1
    assert routed[0].lag_s == 57600.0
    assert routed[0].via_transfer_model == "blood_transit"


def test_topological_order():
    resolver = CausalResolver(_EDGES, transfer_model_ids=[])
    order = resolver.topological_order()
    bm_idx  = order.index("bm_haematopoiesis")
    thy_idx = order.index("thymus_selection")
    pln_idx = order.index("peripheral_ln")
    assert bm_idx < thy_idx < pln_idx
