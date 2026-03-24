"""
Component 3 — Causal Resolver

Builds a directed acyclic graph (DAG) from the config edges list.
Routes export signals from source models to target models, computing
scheduled delivery times from lag specifications.

Lag formats:
  constant:<N>   — fixed N-second delay
  model:<ID>     — lag is read from the transfer model's lag_s field
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

import networkx as nx

logger = logging.getLogger(__name__)

_CONSTANT_RE = re.compile(r"^constant:(\d+)$")
_MODEL_RE    = re.compile(r"^model:([a-zA-Z_][a-zA-Z0-9_]*)$")


@dataclass
class RoutedSignal:
    """A signal queued for delivery to a target model."""
    signal_id: str
    source_model: str
    target_model: str
    payload: dict           # the export_signal dict from the source ISSL
    emitted_at_s: float
    deliver_at_s: float     # scheduled delivery time
    lag_s: float
    via_transfer_model: str | None = None


class CausalResolver:
    """DAG-based signal router."""

    def __init__(self, edges: list[dict], transfer_model_ids: list[str]) -> None:
        self._edges = edges
        self._transfer_model_ids = set(transfer_model_ids or [])
        self._dag = nx.DiGraph()
        self._build_dag(edges)

    # ------------------------------------------------------------------

    def _build_dag(self, edges: list[dict]) -> None:
        for edge in edges:
            src = edge["source"]
            tgt = edge["target"]
            self._dag.add_edge(src, tgt, **edge)

        if not nx.is_directed_acyclic_graph(self._dag):
            raise ValueError(
                "Configuration graph contains a cycle — signal routing would deadlock. "
                f"Edges: {edges}"
            )
        logger.info("Causal DAG built: %d nodes, %d edges",
                    self._dag.number_of_nodes(), self._dag.number_of_edges())

    def get_edges_from(self, model_id: str) -> list[dict]:
        """Return all edge dicts originating from *model_id*."""
        return [
            data
            for _, _, data in self._dag.out_edges(model_id, data=True)
        ]

    def route(
        self,
        model_id: str,
        issl_record: dict,
        sim_time_s: float,
        transfer_issl_records: dict[str, dict],
    ) -> list[RoutedSignal]:
        """
        For each export_signal in *issl_record*, find matching edges and
        create RoutedSignal objects with delivery times resolved.

        *transfer_issl_records* maps transfer_model_id → ISSL record
        (already computed by TransferDispatcher for model edges).
        """
        routed: list[RoutedSignal] = []
        export_signals = issl_record.get("export_signals", [])
        edges = self.get_edges_from(model_id)

        for edge in edges:
            signal_id = edge["signal_id"]
            target = edge["target"]
            lag_spec = edge["lag"]
            threshold = edge.get("activation_threshold") or 0.0

            # Find the matching export signal
            payload = next(
                (s for s in export_signals if s["signal_id"] == signal_id),
                None,
            )
            if payload is None:
                logger.debug("No export signal '%s' in %s ISSL — skipping edge",
                             signal_id, model_id)
                continue

            flux = payload.get("flux", 0.0)
            if flux <= threshold:
                logger.debug("Signal '%s' flux %.4f below threshold %.4f — not routed",
                             signal_id, flux, threshold)
                continue

            # Resolve lag
            lag_s, via_transfer = self._resolve_lag(lag_spec, transfer_issl_records)

            routed.append(RoutedSignal(
                signal_id=signal_id,
                source_model=model_id,
                target_model=target,
                payload=payload,
                emitted_at_s=sim_time_s,
                deliver_at_s=sim_time_s + lag_s,
                lag_s=lag_s,
                via_transfer_model=via_transfer,
            ))
            logger.info("Routed %s → %s (lag=%.0f s, via=%s)",
                        model_id, target, lag_s, via_transfer)

        return routed

    def _resolve_lag(
        self, lag_spec: str, transfer_records: dict[str, dict]
    ) -> tuple[float, str | None]:
        m_const = _CONSTANT_RE.match(lag_spec)
        if m_const:
            return float(m_const.group(1)), None

        m_model = _MODEL_RE.match(lag_spec)
        if m_model:
            transfer_id = m_model.group(1)
            record = transfer_records.get(transfer_id)
            if record is None:
                raise RuntimeError(
                    f"Transfer model '{transfer_id}' has no ISSL record — "
                    "was TransferDispatcher.dispatch() called first?"
                )
            # lag_s is in export_signals[0]
            lag_s = record["export_signals"][0].get("lag_s") or 0.0
            return float(lag_s), transfer_id

        raise ValueError(f"Unrecognised lag spec: {lag_spec!r}")

    def topological_order(self) -> list[str]:
        """Return model IDs in topological order (sources first)."""
        return list(nx.topological_sort(self._dag))
