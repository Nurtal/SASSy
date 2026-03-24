# Adding a New Model to OISA

This guide walks through adding a hypothetical `spleen_marginal_zone` ODE model as a worked example.

## 1. Create the model directory

```
models/spleen_marginal_zone/
├── __init__.py
├── model.py
├── parameters.yaml
└── tests/
    ├── __init__.py
    └── test_spleen.py
```

## 2. Implement `model.py`

Subclass `ModelBase` and implement two methods:

```python
from models._base.model_base import ModelBase

class SpleenMarginalZone(ModelBase):

    MODEL_ID      = "spleen_marginal_zone"
    MODEL_VERSION = "1"
    DELTA_T_S     = 43_200.0  # 12 h

    def __init__(self, port: str, output_dir: Path) -> None:
        super().__init__(
            model_id=self.MODEL_ID,
            model_version=self.MODEL_VERSION,
            formalism="ODE",
            delta_t_s=self.DELTA_T_S,
            port=port,
            output_dir=output_dir,
        )
        # load parameters, initialise state …

    def _step(self, sim_time_s: float, signals: list[dict]) -> None:
        # Integrate ODE by one DELTA_T_S.
        # Parse any incoming signals from signals[].
        …

    def emit_issl(self, sim_time_s: float) -> dict:
        # Return a dict conforming to schemas/issl_v1.schema.json.
        # Use self._make_watchdog() for the watchdog section.
        return { "envelope": …, "continuous_state": …, … }
```

Key rules:
- `_step()` must advance state by exactly `DELTA_T_S` seconds.
- `emit_issl()` must return a dict that validates against `schemas/issl_v1.schema.json`.
- Any flux to route downstream goes in `export_signals` with a stable `signal_id`.
- Transfer models must populate `lag_s` in `export_signals[0]`.

## 3. Write `parameters.yaml`

```yaml
model_id: spleen_marginal_zone
model_version: "1"

initial_conditions:
  MZ_B: 500000   # marginal zone B cells

parameters:
  r_proliferation:
    value: 0.01
    unit: day^-1
    ci_95: [0.008, 0.012]
    source: "Author et al. YYYY, Journal"
```

List every parameter with value, unit, CI, and a literature source.  These values are automatically included in `internal_parameters` of every ISSL record and in `supplementary/table_S1_parameters.csv`.

## 4. Register in a config graph

Add the model to the relevant `configs/run_*.yaml`:

```yaml
models:
  - id: spleen_marginal_zone
    formalism: ODE
    executable: models/spleen_marginal_zone/model.py
    issl_port: "tcp://localhost:5014"
    delta_t_s: 43200

edges:
  - source: thymus_selection
    signal_id: thymus_selection.naive_T_export
    target: spleen_marginal_zone
    lag: "constant:86400"
    activation_threshold: 0.0
```

## 5. Write tests

Three minimum tests (see `models/bm_haematopoiesis/tests/test_bm.py` as a template):

1. `test_step_no_crash` — step with no signals keeps state non-negative.
2. `test_emit_issl_validates` — emitted record validates against `issl_v1.schema.json`.
3. One biological assertion (e.g. pool grows from sub-set-point initial condition).

## 6. Update `supplementary/table_S1_parameters.csv`

Add one row per parameter following the existing format.
