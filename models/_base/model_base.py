"""
Abstract base class for all OISA models.

Each concrete model must implement:
  _step(sim_time_s, signals)  — advance state by one delta_t_s
  emit_issl(sim_time_s)       — return ISSL record dict

The base class owns:
  - ZeroMQ REP server loop (run())
  - ISSL schema validation before every emit
  - Checkpoint file writing via issl_writer
"""

from __future__ import annotations

import json
import logging
from abc import ABC, abstractmethod
from pathlib import Path

import jsonschema
import zmq

from models._base.issl_writer import write_checkpoint

logger = logging.getLogger(__name__)

# Resolve schemas/ directory relative to this file (two levels up from _base/)
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_ISSL_SCHEMA_PATH = _REPO_ROOT / "schemas" / "issl_v1.schema.json"


def _load_issl_schema() -> dict:
    if not _ISSL_SCHEMA_PATH.exists():
        raise FileNotFoundError(f"ISSL schema not found at {_ISSL_SCHEMA_PATH}")
    return json.loads(_ISSL_SCHEMA_PATH.read_text())


class ModelBase(ABC):
    """Base class for OISA model service processes."""

    def __init__(
        self,
        model_id: str,
        model_version: str,
        formalism: str,
        delta_t_s: float | None,
        port: str,
        output_dir: Path,
    ) -> None:
        self.model_id = model_id
        self.model_version = model_version
        self.formalism = formalism
        self.delta_t_s = delta_t_s
        self.output_dir = Path(output_dir)

        self._issl_schema = _load_issl_schema()
        self._validator = jsonschema.Draft202012Validator(self._issl_schema)

        self._ctx = zmq.Context()
        self._sock: zmq.Socket = self._ctx.socket(zmq.REP)
        self._sock.bind(port)

        logger.info("%s bound to %s", model_id, port)

    # ------------------------------------------------------------------
    # Abstract interface

    @abstractmethod
    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        """Advance the model by one delta_t_s.

        *signals* is a list of export-signal dicts from upstream models,
        as routed by the orchestrator.
        """

    @abstractmethod
    def emit_issl(self, sim_time_s: float) -> dict:
        """Return the current state as a validated ISSL record dict."""

    # ------------------------------------------------------------------
    # Server loop

    def run(self) -> None:
        """Block and serve step/emit/shutdown commands from the orchestrator."""
        logger.info("%s entering serve loop", self.model_id)
        while True:
            try:
                msg: dict = self._sock.recv_json()
            except zmq.ZMQError as exc:
                logger.error("ZMQ receive error: %s", exc)
                break

            cmd = msg.get("cmd")

            if cmd == "step":
                sim_time_s: float = msg["sim_time_s"]
                signals: list[dict] = msg.get("signals", [])
                try:
                    self._step(sim_time_s, signals)
                    self._sock.send_json({"status": "ok"})
                except Exception as exc:
                    logger.exception("Error in _step at t=%s", sim_time_s)
                    self._sock.send_json({"status": "error", "message": str(exc)})

            elif cmd == "emit":
                sim_time_s = msg["sim_time_s"]
                try:
                    record = self.emit_issl(sim_time_s)
                    self._validate(record)
                    write_checkpoint(record, self.output_dir / self.model_id, sim_time_s)
                    self._sock.send_json(record)
                except jsonschema.ValidationError as exc:
                    logger.error("ISSL validation failed: %s", exc.message)
                    self._sock.send_json({"status": "error", "message": exc.message})
                except Exception as exc:
                    logger.exception("Error in emit_issl at t=%s", sim_time_s)
                    self._sock.send_json({"status": "error", "message": str(exc)})

            elif cmd == "shutdown":
                logger.info("%s shutting down", self.model_id)
                self._sock.send_json({"status": "bye"})
                break

            else:
                logger.warning("Unknown command: %s", cmd)
                self._sock.send_json({"status": "error", "message": f"unknown cmd: {cmd}"})

        self._sock.close()
        self._ctx.term()

    # ------------------------------------------------------------------
    # Helpers

    def _validate(self, record: dict) -> None:
        errors = list(self._validator.iter_errors(record))
        if errors:
            # Report the first (most informative) error
            raise jsonschema.ValidationError(
                f"ISSL validation failed for {self.model_id}: "
                + "; ".join(e.message for e in errors[:3])
            )

    def _make_watchdog(
        self,
        sim_time_s: float,
        health_status: str = "ok",
        ood_flag: bool = False,
        divergence_score: float | None = None,
    ) -> dict:
        """Build a standard watchdog section."""
        return {
            "health_status": health_status,
            "ood_flag": ood_flag,
            "divergence_score": divergence_score,
            "next_checkpoint_s": sim_time_s + (self.delta_t_s or 0),
        }
