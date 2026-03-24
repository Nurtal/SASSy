"""
Component 1 — ISSL Ingestion

Validates incoming ISSL records against the JSON schema and writes
checkpoint files to disk via issl_writer.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import jsonschema

from models._base.issl_writer import write_checkpoint

logger = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_ISSL_SCHEMA_PATH = _REPO_ROOT / "schemas" / "issl_v1.schema.json"


class ISSLIngestion:
    """Validate and persist ISSL records received from model processes."""

    def __init__(self, output_dir: Path) -> None:
        self._output_dir = Path(output_dir)
        schema_text = _ISSL_SCHEMA_PATH.read_text()
        self._schema = json.loads(schema_text)
        self._validator = jsonschema.Draft202012Validator(self._schema)

    def validate(self, record: dict) -> None:
        """Raise jsonschema.ValidationError if *record* does not conform."""
        errors = list(self._validator.iter_errors(record))
        if errors:
            best = jsonschema.exceptions.best_match(errors)
            raise jsonschema.ValidationError(best.message)

    def ingest(self, model_id: str, record: dict, sim_time_s: float) -> Path:
        """Validate *record* and write it to disk. Returns the path written."""
        self.validate(record)
        model_dir = self._output_dir / "issl" / model_id
        path = write_checkpoint(record, model_dir, sim_time_s)
        logger.debug("Ingested ISSL for %s at t=%s → %s", model_id, sim_time_s, path)
        return path
