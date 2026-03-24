"""
Component 6 — Transfer Dispatcher

When a causal edge specifies lag: "model:<ID>", the orchestrator invokes
the named transfer model via ZMQ before routing the signal.  This component
manages those synchronous step→emit call sequences.
"""

from __future__ import annotations

import logging

import zmq

logger = logging.getLogger(__name__)


class TransferDispatcher:
    """Invoke transfer models on-demand to compute lag_s and transformed signals."""

    def __init__(self, model_sockets: dict[str, zmq.Socket]) -> None:
        """
        *model_sockets* maps model_id → connected ZMQ REQ socket.
        Only transfer model IDs need to be present.
        """
        self._sockets = model_sockets

    def dispatch(
        self,
        transfer_model_id: str,
        source_signal: dict,
        sim_time_s: float,
    ) -> dict:
        """
        Send step + emit commands to the transfer model.

        Returns the ISSL record emitted by the transfer model.  The caller
        (CausalResolver) reads export_signals[0].lag_s from this record.
        """
        sock = self._sockets.get(transfer_model_id)
        if sock is None:
            raise RuntimeError(
                f"No ZMQ socket for transfer model '{transfer_model_id}'. "
                "Check configuration and that the model process is running."
            )

        # Step command: pass the source signal as initial condition
        step_cmd = {
            "cmd":        "step",
            "sim_time_s": sim_time_s,
            "signals":    [source_signal],
        }
        sock.send_json(step_cmd)
        ack = sock.recv_json()
        if ack.get("status") != "ok":
            raise RuntimeError(
                f"Transfer model '{transfer_model_id}' step failed: {ack}"
            )

        # Emit command: retrieve the ISSL record
        emit_cmd = {"cmd": "emit", "sim_time_s": sim_time_s}
        sock.send_json(emit_cmd)
        record = sock.recv_json()

        if "status" in record and record["status"] == "error":
            raise RuntimeError(
                f"Transfer model '{transfer_model_id}' emit failed: {record}"
            )

        logger.info(
            "TransferDispatcher: '%s' → lag_s=%.1f, delivered=%.2f",
            transfer_model_id,
            record["export_signals"][0].get("lag_s", 0.0),
            record["export_signals"][0].get("flux", 0.0),
        )
        return record
