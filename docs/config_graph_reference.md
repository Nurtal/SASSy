# Configuration Graph Field Reference — v1

The configuration graph is a YAML file validated against `schemas/config_graph_v1.schema.json`.

## Top-level fields

| Field | Required | Description |
|---|---|---|
| `oisa_version` | yes | Schema version string (e.g. `"1.0"`) |
| `models` | yes | List of model service declarations |
| `edges` | no | Signal routing edges between models |
| `transfer_models` | no | Model IDs that act as transfer models (on-demand) |
| `global_clock` | yes | Simulation time bounds and OISSL checkpoint interval |
| `renderer` | no | Optional 3D render stream configuration |

---

## `models[]`

| Field | Required | Description |
|---|---|---|
| `id` | yes | Unique identifier (alphanumeric + underscore) |
| `formalism` | yes | `"ODE"` \| `"ABM"` \| `"hybrid"` |
| `executable` | yes | Path to the model's Python entry point |
| `issl_port` | yes | ZMQ bind/connect address (`tcp://localhost:<port>`) |
| `delta_t_s` | no | Checkpoint interval in seconds; `null` for on-demand transfer models |

---

## `edges[]`

| Field | Required | Description |
|---|---|---|
| `source` | yes | Sending model `id` |
| `signal_id` | yes | Must match a `signal_id` in the source model's `export_signals` |
| `target` | yes | Receiving model `id` |
| `lag` | yes | Delay specification (see below) |
| `activation_threshold` | no | Minimum flux to trigger routing (default 0) |

### `lag` formats

| Pattern | Example | Description |
|---|---|---|
| `constant:<N>` | `constant:172800` | Fixed N-second delay |
| `model:<ID>` | `model:blood_transit` | Lag computed by the named transfer model via its `export_signals[0].lag_s` |

---

## `global_clock`

| Field | Required | Description |
|---|---|---|
| `start_s` | yes | Simulation start time in seconds |
| `end_s` | yes | Simulation end time in seconds |
| `checkpoint_interval_s` | yes | Interval at which OISSL records are emitted |

Note: per-model `delta_t_s` and `checkpoint_interval_s` are independent.  A model with `delta_t_s: 21600` fires 4× per OISSL checkpoint of `86400`.

---

## `renderer`

| Field | Description |
|---|---|
| `format` | Output format tag (e.g. `"OISA-render-v1"`) |
| `target` | `file:<path>` or `tcp://<host>:<port>` |
| `emit_interval_s` | How often to write to the render stream |
| `scene_schema_uri` | Path to Godot scene schema JSON |
