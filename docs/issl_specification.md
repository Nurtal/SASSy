# ISSL Field Reference — v1

Internal Simulation State Log (ISSL) is the JSON-LD record emitted by each OISA model at every checkpoint.  Schema: `schemas/issl_v1.schema.json`.

## Top-level structure

```json
{
  "envelope":             { ... },
  "continuous_state":     [ ... ],
  "discrete_events":      [ ... ],
  "export_signals":       [ ... ],
  "internal_parameters":  [ ... ],
  "watchdog":             { ... }
}
```

---

## `envelope`

| Field | Type | Description |
|---|---|---|
| `model_id` | string | Stable identifier matching the config graph |
| `model_version` | string | Semantic version of the model implementation |
| `sim_time_s` | number ≥ 0 | GSimT at the *end* of the step that produced this record |
| `schema_uri` | string | URI reference for the schema (e.g. `schemas/issl_v1.schema.json`) |
| `formalism` | `"ODE"` \| `"ABM"` \| `"hybrid"` | Computational formalism |

---

## `continuous_state`

Array of entity snapshots.  Each entry:

| Field | Type | Description |
|---|---|---|
| `entity_id` | string | OBO URI (e.g. `CL:0000037`) or `MODEL_X_custom:name` |
| `entity_class` | `"obo"` \| `"custom"` | Whether the ID maps to an OBO term |
| `label` | string | Human-readable name |
| `count` | number | Cell count (or concentration) at sim_time_s |
| `unit` | string | e.g. `"cells"` or `"cells·mL^-1"` |
| `fitness` | number \| null | Scalar fitness (ABM) or null |
| `surface_markers` | string[] | CD marker profile |
| `ci_95` | [lo, hi] | 95 % credible interval; ODE: analytical sensitivity; ABM: empirical across realisations |

---

## `discrete_events`

Array of punctual events that occurred during the step:

| Field | Type | Description |
|---|---|---|
| `event_type` | string | e.g. `"positive_selection"`, `"apoptosis"` |
| `entity_id` | string | OBO term for the process (e.g. `GO:0045058`) |
| `count` | integer ≥ 0 | Number of events |
| `sim_time_s` | number | GSimT when the events occurred |
| `n_realisations` | integer | Number of ABM realisations contributing (ODE: omit) |
| `variance` | number \| null | Empirical variance across realisations |

---

## `export_signals`

Fluxes available for routing to downstream models:

| Field | Type | Description |
|---|---|---|
| `signal_id` | string | Matches the `signal_id` on a config graph edge |
| `entity_id` | string | Entity being exported |
| `flux` | number | Export rate in `unit` |
| `unit` | string | e.g. `"cells·day^-1"` or `"cells·checkpoint^-1"` |
| `lag_s` | number \| null | For transfer models only: computed transit time in seconds |
| `ci_95` | [lo, hi] | 95 % CI on the flux value |

> **Transfer models** must populate `lag_s`.  The orchestrator reads this field via `TransferDispatcher` to schedule signal delivery.

---

## `internal_parameters`

All kinetic parameters used during the step:

| Field | Type | Description |
|---|---|---|
| `param_id` | string | Key matching `parameters.yaml` |
| `value` | number | Current value |
| `unit` | string | |
| `ci_95` | [lo, hi] | Posterior credible interval |
| `source` | string | Literature reference |

---

## `watchdog`

Health status at the end of the step:

| Field | Type | Description |
|---|---|---|
| `health_status` | `"ok"` \| `"degraded"` \| `"failed"` | Overall health |
| `ood_flag` | boolean | True if any state variable is outside expected biological range |
| `divergence_score` | number \| null | Scalar measure of deviation from expected steady state |
| `next_checkpoint_s` | number | GSimT of this model's next step |

The orchestrator's `Watchdog` component fires a `watchdog_alert` composition event when `health_status != "ok"` or `ood_flag == true`.
