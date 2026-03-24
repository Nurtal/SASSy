"""
Component 2 — Scheduler (Global Simulation Time — GSimT)

Maintains a priority queue of per-model checkpoint events.  Models with
different delta_t_s fire at different rates; the global OISSL checkpoint
fires at checkpoint_interval_s regardless.

The scheduler drives the main simulation loop by:
  1. Popping the next event (earliest sim_time_s)
  2. Sending step + emit commands to the model via ZMQ
  3. Collecting the ISSL record and passing it to the pipeline
  4. Enqueueing the model's next checkpoint
  5. Emitting an OISSL record when the global clock ticks

Pending signals (routed by CausalResolver) are delivered when
deliver_at_s <= current_model_time for the target model.
"""

from __future__ import annotations

import heapq
import logging
from dataclasses import dataclass, field
from typing import Callable

import zmq

from orchestrator.components.causal_resolver import CausalResolver, RoutedSignal
from orchestrator.components.constraint_engine import ConstraintEngine
from orchestrator.components.ingestion import ISSLIngestion
from orchestrator.components.output_aggregator import OutputAggregator
from orchestrator.components.state_registry import StateRegistry
from orchestrator.components.transfer_dispatcher import TransferDispatcher
from orchestrator.components.watchdog import Watchdog

logger = logging.getLogger(__name__)


@dataclass(order=True)
class _Event:
    sim_time_s: float
    model_id: str = field(compare=False)
    is_global_checkpoint: bool = field(compare=False, default=False)


class Scheduler:
    """Global simulation clock and event-driven loop."""

    def __init__(
        self,
        global_clock: dict,
        model_cfgs: list[dict],
        model_sockets: dict[str, zmq.Socket],
        ingestion: ISSLIngestion,
        causal_resolver: CausalResolver,
        constraint_engine: ConstraintEngine,
        state_registry: StateRegistry,
        transfer_dispatcher: TransferDispatcher,
        output_aggregator: OutputAggregator,
        watchdog: Watchdog,
    ) -> None:
        self._start_s: float = float(global_clock["start_s"])
        self._end_s: float   = float(global_clock["end_s"])
        self._ckpt_interval: float = float(global_clock["checkpoint_interval_s"])

        self._model_cfgs: dict[str, dict] = {m["id"]: m for m in model_cfgs}
        self._sockets = model_sockets

        self._ingestion   = ingestion
        self._resolver    = causal_resolver
        self._constraints = constraint_engine
        self._registry    = state_registry
        self._dispatcher  = transfer_dispatcher
        self._aggregator  = output_aggregator
        self._watchdog    = watchdog

        # Current GSimT
        self._current_s: float = self._start_s

        # Pending routed signals waiting for delivery
        self._pending_signals: list[RoutedSignal] = []

        # Composition events accumulated since last OISSL checkpoint
        self._composition_events: list[dict] = []

        # Priority queue
        self._queue: list[_Event] = []
        self._build_initial_queue()

    # ------------------------------------------------------------------

    def _build_initial_queue(self) -> None:
        for model_id, cfg in self._model_cfgs.items():
            delta_t = cfg.get("delta_t_s")
            if delta_t is None:
                continue  # transfer models are invoked on-demand
            heapq.heappush(self._queue, _Event(
                sim_time_s=self._start_s + delta_t,
                model_id=model_id,
            ))

        # First global OISSL checkpoint
        heapq.heappush(self._queue, _Event(
            sim_time_s=self._start_s + self._ckpt_interval,
            model_id="__global__",
            is_global_checkpoint=True,
        ))
        logger.info("Scheduler: queue initialised with %d events", len(self._queue))

    # ------------------------------------------------------------------

    def run(self) -> None:
        """Main simulation loop. Blocks until end_s is reached."""
        logger.info("Simulation start: t=%.0f s → t=%.0f s (%.1f days)",
                    self._start_s, self._end_s,
                    (self._end_s - self._start_s) / 86_400)

        while self._queue:
            event = heapq.heappop(self._queue)
            if event.sim_time_s > self._end_s:
                break

            self._current_s = event.sim_time_s

            if event.is_global_checkpoint:
                self._emit_oissl(event.sim_time_s)
                next_t = event.sim_time_s + self._ckpt_interval
                if next_t <= self._end_s:
                    heapq.heappush(self._queue, _Event(
                        sim_time_s=next_t,
                        model_id="__global__",
                        is_global_checkpoint=True,
                    ))
            else:
                self._step_model(event.model_id, event.sim_time_s)

        logger.info("Simulation complete at GSimT=%.0f s", self._current_s)

    # ------------------------------------------------------------------

    def _step_model(self, model_id: str, sim_time_s: float) -> None:
        """Send step + emit to *model_id*, process the returned ISSL."""
        cfg = self._model_cfgs[model_id]
        delta_t = cfg["delta_t_s"]
        sock = self._sockets[model_id]

        # Collect pending signals whose delivery time has arrived
        due_signals = [
            rs for rs in self._pending_signals
            if rs.target_model == model_id and rs.deliver_at_s <= sim_time_s
        ]
        self._pending_signals = [
            rs for rs in self._pending_signals
            if not (rs.target_model == model_id and rs.deliver_at_s <= sim_time_s)
        ]
        signal_payloads = [rs.payload for rs in due_signals]

        # ---- step ----
        step_t = sim_time_s - delta_t  # step starts one delta_t ago
        sock.send_json({"cmd": "step", "sim_time_s": step_t, "signals": signal_payloads})
        ack = sock.recv_json()
        if ack.get("status") != "ok":
            logger.error("Model %s step error at t=%.0f: %s", model_id, sim_time_s, ack)

        # ---- emit ----
        sock.send_json({"cmd": "emit", "sim_time_s": sim_time_s})
        record = sock.recv_json()

        if "status" in record and record["status"] == "error":
            logger.error("Model %s emit error at t=%.0f: %s", model_id, sim_time_s, record)
        else:
            self._process_issl(model_id, record, sim_time_s)

        # Enqueue next checkpoint for this model
        next_t = sim_time_s + delta_t
        if next_t <= self._end_s:
            heapq.heappush(self._queue, _Event(sim_time_s=next_t, model_id=model_id))

    # ------------------------------------------------------------------

    def _process_issl(self, model_id: str, record: dict, sim_time_s: float) -> None:
        """Full ISSL processing pipeline for one received record."""
        # 1. Validate & persist
        try:
            self._ingestion.ingest(model_id, record, sim_time_s)
        except Exception as exc:
            logger.error("Ingestion failed for %s at t=%.0f: %s", model_id, sim_time_s, exc)

        # 2. Watchdog
        alert = self._watchdog.check(model_id, record, sim_time_s)
        if alert:
            events = self._watchdog.to_composition_events([alert])
            self._composition_events.extend(events)

        # 3. Constraint checks
        violations = self._constraints.check(model_id, record)
        if violations:
            events = self._constraints.to_composition_events(violations, sim_time_s)
            self._composition_events.extend(events)

        # 4. State registry update
        self._registry.update(model_id, record)

        # 5. Signal routing via causal resolver
        transfer_records: dict[str, dict] = {}
        for edge in self._resolver.get_edges_from(model_id):
            lag_spec = edge.get("lag", "constant:0")
            if lag_spec.startswith("model:"):
                transfer_id = lag_spec[6:]
                export_signals = record.get("export_signals", [])
                if export_signals:
                    try:
                        t_record = self._dispatcher.dispatch(
                            transfer_id, export_signals[0], sim_time_s
                        )
                        transfer_records[transfer_id] = t_record
                        self._composition_events.append({
                            "type":         "transfer_invoked",
                            "sim_time_s":   sim_time_s,
                            "source_model": model_id,
                            "target_model": edge["target"],
                            "signal_id":    edge["signal_id"],
                            "lag_s":        t_record["export_signals"][0].get("lag_s"),
                            "value":        t_record["export_signals"][0].get("flux"),
                        })
                    except Exception as exc:
                        logger.error("TransferDispatcher error: %s", exc)

        routed = self._resolver.route(model_id, record, sim_time_s, transfer_records)
        for rs in routed:
            self._pending_signals.append(rs)
            self._composition_events.append({
                "type":         "signal_emitted",
                "sim_time_s":   sim_time_s,
                "signal_id":    rs.signal_id,
                "source_model": rs.source_model,
                "target_model": rs.target_model,
                "lag_s":        rs.lag_s,
                "value":        rs.payload.get("flux"),
            })

        # Record deliveries that just happened (due_signals from _step_model)
        for rs in [
            rs for rs in self._pending_signals
            if rs.target_model == model_id and rs.deliver_at_s <= sim_time_s
        ]:
            self._composition_events.append({
                "type":         "signal_delivered",
                "sim_time_s":   sim_time_s,
                "signal_id":    rs.signal_id,
                "source_model": rs.source_model,
                "target_model": rs.target_model,
                "lag_s":        rs.lag_s,
                "value":        rs.payload.get("flux"),
            })

    # ------------------------------------------------------------------

    def _emit_oissl(self, sim_time_s: float) -> None:
        """Assemble and write an OISSL record at the global checkpoint."""
        self._aggregator.emit(
            sim_time_s=sim_time_s,
            gis_snapshot=self._registry.snapshot(),
            composition_events=list(self._composition_events),
            render_entities=self._registry.render_entities(),
        )
        self._composition_events.clear()
        logger.info("OISSL emitted at GSimT=%.0f s (day %.1f)",
                    sim_time_s, sim_time_s / 86_400)

    @property
    def current_time(self) -> float:
        return self._current_s
