"""
Component 5 — State Registry

Maintains the Global Immune State (GIS): a merged snapshot of all model
continuous_state arrays, keyed by entity_id.  The latest value from each
model is kept; if two models track the same entity_id the most recently
updated wins.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class EntitySnapshot:
    entity_id: str
    label: str
    count: float
    unit: str
    model_id: str
    sim_time_s: float
    ci_95: list[float] | None = None
    surface_markers: list[str] = field(default_factory=list)


class StateRegistry:
    """Global Immune State registry."""

    def __init__(self) -> None:
        # entity_id → EntitySnapshot
        self._state: dict[str, EntitySnapshot] = {}
        # model_id → list of entity snapshots (for OISSL per-model sections)
        self._model_snapshots: dict[str, list[dict]] = {}

    def update(self, model_id: str, record: dict) -> None:
        """Merge *record*'s continuous_state into the GIS."""
        sim_time_s = record["envelope"]["sim_time_s"]
        snapshot_list = record.get("continuous_state", [])

        self._model_snapshots[model_id] = snapshot_list

        for entity in snapshot_list:
            eid = entity["entity_id"]
            self._state[eid] = EntitySnapshot(
                entity_id=eid,
                label=entity.get("label", ""),
                count=entity["count"],
                unit=entity.get("unit", "cells"),
                model_id=model_id,
                sim_time_s=sim_time_s,
                ci_95=entity.get("ci_95"),
                surface_markers=entity.get("surface_markers", []),
            )

        logger.debug("StateRegistry updated for %s (%d entities)", model_id, len(snapshot_list))

    def snapshot(self) -> list[dict]:
        """Return the GIS as a list of per-model dicts for OISSL global_immune_state."""
        result = []
        for model_id, entities in self._model_snapshots.items():
            sim_time_s = max(
                (e.get("sim_time_s", 0.0) for e in entities),
                default=0.0,
            ) if entities else 0.0
            # Use the most recent sim_time_s from the registry for this model
            for snap in self._state.values():
                if snap.model_id == model_id:
                    sim_time_s = snap.sim_time_s
                    break
            result.append({
                "model_id":   model_id,
                "sim_time_s": sim_time_s,
                "snapshot":   entities,
            })
        return result

    def get(self, entity_id: str) -> EntitySnapshot | None:
        return self._state.get(entity_id)

    def render_entities(self) -> list[dict]:
        """Return render_entities list for OISSL (one entry per tracked entity)."""
        return [
            {
                "entity_id": snap.entity_id,
                "label":     snap.label,
                "model_id":  snap.model_id,
                "count":     snap.count,
            }
            for snap in self._state.values()
        ]
