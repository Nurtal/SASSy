"""
Component 8 — Output Aggregator

Assembles OISSL records from the StateRegistry, composition_events list,
and render_entities, then writes them to disk (and optionally to a TCP
socket for the Godot renderer).
"""

from __future__ import annotations

import json
import logging
import socket
import time
from pathlib import Path

import jsonschema

logger = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_OISSL_SCHEMA_PATH = _REPO_ROOT / "schemas" / "oissl_v1.schema.json"

ORCHESTRATOR_VERSION = "0.1.0"


class OutputAggregator:
    """Assemble and persist OISSL records."""

    def __init__(
        self,
        output_dir: Path,
        run_id: str,
        config_uri: str,
        renderer_cfg: dict | None = None,
    ) -> None:
        self._output_dir = Path(output_dir)
        self._run_id = run_id
        self._config_uri = config_uri
        self._renderer_cfg = renderer_cfg or {}

        schema_text = _OISSL_SCHEMA_PATH.read_text()
        self._schema = json.loads(schema_text)
        self._validator = jsonschema.Draft202012Validator(self._schema)

        self._oissl_dir = self._output_dir / "oissl"
        self._oissl_dir.mkdir(parents=True, exist_ok=True)

        # Optional render stream (file or TCP)
        self._render_file: Path | None = None
        self._render_sock: socket.socket | None = None
        self._setup_renderer()

    def _setup_renderer(self) -> None:
        target = self._renderer_cfg.get("target", "")
        if target.startswith("file:"):
            self._render_file = Path(target[5:])
            self._render_file.parent.mkdir(parents=True, exist_ok=True)
        elif target.startswith("tcp://"):
            # Best-effort TCP connection to renderer
            try:
                host, port_str = target[6:].rsplit(":", 1)
                self._render_sock = socket.create_connection((host, int(port_str)), timeout=2.0)
                logger.info("OutputAggregator: render stream connected to %s", target)
            except OSError as exc:
                logger.warning("OutputAggregator: renderer not reachable (%s) — skipping", exc)
                self._render_sock = None

    def emit(
        self,
        sim_time_s: float,
        gis_snapshot: list[dict],
        composition_events: list[dict],
        render_entities: list[dict],
    ) -> Path:
        """Assemble, validate, and write one OISSL record. Returns the path written."""
        record = {
            "envelope": {
                "orchestrator_version": ORCHESTRATOR_VERSION,
                "sim_time_s":           sim_time_s,
                "run_id":               self._run_id,
                "config_uri":           self._config_uri,
            },
            "global_immune_state": gis_snapshot,
            "composition_events":  composition_events,
            "render_entities":     render_entities,
        }

        # Validate
        errors = list(self._validator.iter_errors(record))
        if errors:
            best = jsonschema.exceptions.best_match(errors)
            logger.error("OISSL validation error at t=%s: %s", sim_time_s, best.message)
            # Write anyway — don't lose data
        path = self._oissl_dir / f"checkpoint_{int(sim_time_s)}.json"
        path.write_text(json.dumps(record, indent=2))
        logger.info("OISSL checkpoint written: %s", path.name)

        # Render stream
        self._write_render(record)

        return path

    def _write_render(self, record: dict) -> None:
        line = json.dumps(record) + "\n"
        if self._render_file:
            with self._render_file.open("a") as f:
                f.write(line)
        if self._render_sock:
            try:
                self._render_sock.sendall(line.encode())
            except OSError as exc:
                logger.warning("Render stream send failed: %s", exc)
                self._render_sock = None

    def flush(self) -> None:
        """Close render stream connections at end of run."""
        if self._render_sock:
            try:
                self._render_sock.close()
            except OSError:
                pass
