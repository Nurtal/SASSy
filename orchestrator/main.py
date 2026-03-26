#!/usr/bin/env python3
"""
OISA Orchestrator — entry point

Reads a configuration graph (YAML), launches model subprocesses, connects
to them via ZMQ, and drives the simulation loop.

Usage:
  python orchestrator/main.py \
      --config configs/run_COMP3_full_graph.yaml \
      --output-dir logs/COMP3 \
      --log-level INFO
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import jsonschema
import yaml
import zmq

from orchestrator.components.calibration_bridge import CalibrationBridge
from orchestrator.components.causal_resolver import CausalResolver
from orchestrator.components.constraint_engine import ConstraintEngine
from orchestrator.components.ingestion import ISSLIngestion
from orchestrator.components.output_aggregator import OutputAggregator
from orchestrator.components.scheduler import Scheduler
from orchestrator.components.state_registry import StateRegistry
from orchestrator.components.transfer_dispatcher import TransferDispatcher
from orchestrator.components.watchdog import Watchdog

logger = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CONFIG_SCHEMA_PATH = _REPO_ROOT / "schemas" / "config_graph_v1.schema.json"

_MODEL_CONNECT_TIMEOUT_S = 10.0   # seconds to wait for each model to come up
_MODEL_PING_INTERVAL_MS  = 200


def _setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(name)-32s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )


def _load_config(config_path: Path) -> dict:
    """Parse YAML config and validate against config_graph_v1.schema.json."""
    cfg = yaml.safe_load(config_path.read_text())
    schema = json.loads(_CONFIG_SCHEMA_PATH.read_text())
    try:
        jsonschema.validate(cfg, schema)
    except jsonschema.ValidationError as exc:
        logger.error("Configuration validation failed: %s", exc.message)
        sys.exit(1)
    return cfg


def _launch_models(
    model_cfgs: list[dict], output_dir: Path
) -> dict[str, subprocess.Popen]:
    """Start each model as a subprocess. Returns {model_id: Popen}."""
    procs: dict[str, subprocess.Popen] = {}
    for cfg in model_cfgs:
        model_id   = cfg["id"]
        executable = cfg["executable"]
        port       = cfg["issl_port"].replace("localhost", "*")
        model_out  = output_dir / "issl" / model_id
        model_out.mkdir(parents=True, exist_ok=True)

        extra_args = cfg.get("model_args", [])
        cmd = [
            sys.executable, executable,
            "--port", port,
            "--output-dir", str(output_dir / "issl"),
            *extra_args,
        ]
        logger.info("Launching model %s: %s", model_id, " ".join(cmd))
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        procs[model_id] = proc
    return procs


def _connect_sockets(
    model_cfgs: list[dict], ctx: zmq.Context
) -> dict[str, zmq.Socket]:
    """Connect REQ sockets to each model and wait for them to respond."""
    sockets: dict[str, zmq.Socket] = {}

    for cfg in model_cfgs:
        model_id = cfg["id"]
        port     = cfg["issl_port"]

        sock: zmq.Socket = ctx.socket(zmq.REQ)
        sock.setsockopt(zmq.LINGER, 0)
        sock.connect(port)
        sockets[model_id] = sock

    # Wait for all models to come up by sending a benign "ping" (step at t=-1)
    for model_id, port in [(cfg["id"], cfg["issl_port"]) for cfg in model_cfgs]:
        logger.info("Waiting for model %s …", model_id)
        connected = False
        deadline = time.monotonic() + _MODEL_CONNECT_TIMEOUT_S
        while time.monotonic() < deadline:
            # Recreate socket each attempt: REQ sockets cannot send again after
            # a send that was never replied to.
            sock = ctx.socket(zmq.REQ)
            sock.setsockopt(zmq.LINGER, 0)
            sock.connect(port)
            sock.send_json({"cmd": "step", "sim_time_s": -1.0, "signals": []})
            if sock.poll(_MODEL_PING_INTERVAL_MS):
                reply = sock.recv_json()
                sock.close()
                connected = True
                logger.info("Model %s ready (reply: %s)", model_id, reply.get("status"))
                break
            sock.close()
            time.sleep(0.05)
        if not connected:
            logger.error("Model %s did not respond within %.0f s — aborting.",
                         model_id, _MODEL_CONNECT_TIMEOUT_S)
            sys.exit(1)
        # Re-open the final persistent socket for this model
        sock = ctx.socket(zmq.REQ)
        sock.setsockopt(zmq.LINGER, 0)
        sock.connect(port)
        sockets[model_id] = sock

    return sockets


def _shutdown_models(
    procs: dict[str, subprocess.Popen],
    sockets: dict[str, zmq.Socket],
) -> None:
    for model_id, sock in sockets.items():
        try:
            sock.send_json({"cmd": "shutdown"})
            if sock.poll(1000):
                sock.recv_json()
        except zmq.ZMQError:
            pass
        sock.close()
    for model_id, proc in procs.items():
        try:
            proc.wait(timeout=5.0)
        except subprocess.TimeoutExpired:
            proc.terminate()
            proc.wait()
        logger.info("Model %s exited with code %s", model_id, proc.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(description="OISA Orchestrator")
    parser.add_argument("--config",     required=True,  help="Path to config graph YAML")
    parser.add_argument("--output-dir", default="logs", help="Base output directory")
    parser.add_argument("--log-level",  default="INFO", help="Logging level")
    parser.add_argument("--run-id",     default=None,   help="Run identifier (default: timestamp)")
    args = parser.parse_args()

    _setup_logging(args.log_level)

    config_path = Path(args.config)
    if not config_path.exists():
        logger.error("Config file not found: %s", config_path)
        sys.exit(1)

    cfg = _load_config(config_path)

    run_id = args.run_id or datetime.now().strftime("%Y%m%dT%H%M%S")
    output_dir = Path(args.output_dir) / run_id
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Run ID: %s  →  %s", run_id, output_dir)

    model_cfgs        = cfg["models"]
    edges             = cfg.get("edges", [])
    transfer_ids      = cfg.get("transfer_models", [])
    global_clock      = cfg["global_clock"]
    renderer_cfg      = cfg.get("renderer", {})

    # --- Launch model subprocesses ---
    procs = _launch_models(model_cfgs, output_dir)

    # --- ZMQ context and sockets ---
    ctx = zmq.Context()
    try:
        sockets = _connect_sockets(model_cfgs, ctx)
    except SystemExit:
        for p in procs.values():
            p.terminate()
        ctx.term()
        raise

    # --- Instantiate components ---
    ingestion   = ISSLIngestion(output_dir=output_dir)
    resolver    = CausalResolver(edges=edges, transfer_model_ids=transfer_ids)
    constraints = ConstraintEngine()
    registry    = StateRegistry()
    dispatcher  = TransferDispatcher(model_sockets={
        tid: sockets[tid] for tid in transfer_ids if tid in sockets
    })
    aggregator  = OutputAggregator(
        output_dir=output_dir,
        run_id=run_id,
        config_uri=str(config_path),
        renderer_cfg=renderer_cfg,
    )
    watchdog    = Watchdog()

    scheduler = Scheduler(
        global_clock=global_clock,
        model_cfgs=model_cfgs,
        model_sockets=sockets,
        ingestion=ingestion,
        causal_resolver=resolver,
        constraint_engine=constraints,
        state_registry=registry,
        transfer_dispatcher=dispatcher,
        output_aggregator=aggregator,
        watchdog=watchdog,
    )

    # --- Run ---
    try:
        scheduler.run()
    except KeyboardInterrupt:
        logger.info("Interrupted by user.")
    finally:
        aggregator.flush()
        watchdog.report()
        _shutdown_models(procs, sockets)
        ctx.term()
        logger.info("Orchestrator shut down cleanly.")


if __name__ == "__main__":
    main()
