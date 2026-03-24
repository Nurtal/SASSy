"""
Utility: write an ISSL record to disk as a checkpoint file.

Naming convention: checkpoint_<sim_time_s_as_int>.json
Used by ModelBase and by the orchestrator's ingestion component.
"""

from __future__ import annotations

import json
from pathlib import Path


def write_checkpoint(record: dict, output_dir: Path, sim_time_s: float) -> Path:
    """Serialise *record* to ``output_dir/checkpoint_<sim_time_s>.json``.

    Creates *output_dir* if it does not exist.  Returns the path written.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"checkpoint_{int(sim_time_s)}.json"
    path.write_text(json.dumps(record, indent=2))
    return path
