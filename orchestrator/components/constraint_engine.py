"""
Component 4 — Constraint Engine

Checks ISSL records for biological plausibility.  Each constraint is a
named rule applied to an entity or ratio in the continuous_state.
Violations are returned as dicts (not exceptions) and recorded in the
OISSL composition_events section.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ConstraintViolation:
    constraint_id: str
    model_id: str
    entity_id: str
    value: float
    message: str
    severity: str  # "warning" | "error"


# Built-in plausibility rules
# Format: (constraint_id, entity_id, min_val, max_val, severity)
_DEFAULT_RULES: list[tuple[str, str, float, float, str]] = [
    ("HSC_pool_bounds",   "CL:0000037",  500.0,    50_000.0, "warning"),
    ("CLP_non_negative",  "CL:0000051",    0.0,   100_000.0, "warning"),
    ("DN1_non_negative",  "CL:0002420",    0.0,    50_000.0, "warning"),
    ("DP_non_negative",   "CL:0000893",    0.0, 5_000_000.0, "warning"),
    ("CD4_pool_bounds",   "CL:0000624",    0.0, 2_000_000.0, "warning"),
    ("CD8_pool_bounds",   "CL:0000625",    0.0, 1_000_000.0, "warning"),
    ("naive_T_export",    "CL:0000898",    0.0,   100_000.0, "warning"),
]

# CD4/CD8 ratio constraint (applied separately)
_CD4_CD8_RATIO_MIN = 0.5
_CD4_CD8_RATIO_MAX = 10.0


class ConstraintEngine:
    """Biological plausibility checker."""

    def __init__(self, extra_rules: list[dict] | None = None) -> None:
        # Build rule lookup: entity_id → list of (constraint_id, min, max, severity)
        self._rules: dict[str, list[tuple[str, float, float, str]]] = {}
        for cid, eid, lo, hi, sev in _DEFAULT_RULES:
            self._rules.setdefault(eid, []).append((cid, lo, hi, sev))

        if extra_rules:
            for rule in extra_rules:
                eid = rule["entity_id"]
                self._rules.setdefault(eid, []).append((
                    rule["constraint_id"],
                    float(rule.get("min", 0.0)),
                    float(rule.get("max", float("inf"))),
                    rule.get("severity", "warning"),
                ))

    def check(self, model_id: str, record: dict) -> list[ConstraintViolation]:
        """Return a list of violations found in *record*."""
        violations: list[ConstraintViolation] = []

        entity_counts: dict[str, float] = {}
        for entity in record.get("continuous_state", []):
            eid = entity["entity_id"]
            count = entity["count"]
            entity_counts[eid] = count

            if eid in self._rules:
                for cid, lo, hi, sev in self._rules[eid]:
                    if not (lo <= count <= hi):
                        msg = (f"{entity.get('label', eid)} count {count:.1f} "
                               f"outside [{lo}, {hi}]")
                        violations.append(ConstraintViolation(
                            constraint_id=cid,
                            model_id=model_id,
                            entity_id=eid,
                            value=count,
                            message=msg,
                            severity=sev,
                        ))
                        logger.warning("[%s] Constraint %s: %s", model_id, cid, msg)

        # CD4/CD8 ratio check
        cd4 = entity_counts.get("CL:0000624")
        cd8 = entity_counts.get("CL:0000625")
        if cd4 is not None and cd8 is not None and cd8 > 0:
            ratio = cd4 / cd8
            if not (_CD4_CD8_RATIO_MIN <= ratio <= _CD4_CD8_RATIO_MAX):
                msg = f"CD4/CD8 ratio {ratio:.2f} outside [{_CD4_CD8_RATIO_MIN}, {_CD4_CD8_RATIO_MAX}]"
                violations.append(ConstraintViolation(
                    constraint_id="cd4_cd8_ratio",
                    model_id=model_id,
                    entity_id="ratio",
                    value=ratio,
                    message=msg,
                    severity="warning",
                ))
                logger.warning("[%s] Constraint cd4_cd8_ratio: %s", model_id, msg)

        return violations

    def to_composition_events(
        self, violations: list[ConstraintViolation], sim_time_s: float
    ) -> list[dict]:
        """Convert violations to OISSL composition_event dicts."""
        return [
            {
                "type":          "constraint_violation",
                "sim_time_s":    sim_time_s,
                "constraint_id": v.constraint_id,
                "source_model":  v.model_id,
                "severity":      v.severity,
                "message":       v.message,
                "value":         v.value,
            }
            for v in violations
        ]
