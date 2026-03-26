# OISA — Orchestrated Immune Simulation Architecture

**Reference implementation for the paper:**
> *Simulation as a Service: A Formalism-Agnostic Orchestration Framework for Modular Immune Disease Modelling*
> npj Digital Medicine, 2026 — [DOI pending]

---

## Table of contents

1. [What this repository is](#what-this-repository-is)
2. [Repository structure](#repository-structure)
3. [Quick start](#quick-start)
4. [Architecture overview](#architecture-overview)
   - [The ISSL format](#the-issl-format)
   - [The configuration graph](#the-configuration-graph)
   - [The orchestrator](#the-orchestrator)
5. [Models](#models)
   - [Model 1 — Bone marrow haematopoiesis (ODE)](#model-1--bone-marrow-haematopoiesis-ode)
   - [Model 2 — Blood transit kinetics (ODE)](#model-2--blood-transit-kinetics-ode)
   - [Model 3 — Thymic T-cell selection (ABM)](#model-3--thymic-t-cell-selection-abm)
   - [Model 4 — Peripheral lymph node dynamics (ODE)](#model-4--peripheral-lymph-node-dynamics-ode)
6. [Running the simulations](#running-the-simulations)
   - [Simulation runs from the paper](#simulation-runs-from-the-paper)
   - [Writing your own configuration graph](#writing-your-own-configuration-graph)
7. [Output logs and analysis](#output-logs-and-analysis)
   - [ISSL per-model logs](#issl-per-model-logs)
   - [OISSL orchestrator logs](#oissl-orchestrator-logs)
   - [Render stream](#render-stream)
8. [Extending OISA](#extending-oisa)
   - [Adding a new model](#adding-a-new-model)
   - [Adding a transfer model on an edge](#adding-a-transfer-model-on-an-edge)
   - [Using custom (non-ontologised) entities](#using-custom-non-ontologised-entities)
9. [Schemas and specifications](#schemas-and-specifications)
10. [Related repositories](#related-repositories)
11. [Contributing](#contributing)
12. [Citation](#citation)
13. [Licence](#licence)

---

## What this repository is

OISA is an architecture for composing independently developed computational models of immune system components — regardless of whether those models are built as ordinary differential equations (ODEs), agent-based models (ABMs), or hybrid approaches — into a single coordinated simulation.

The core idea is simple: **each model is a service**. It runs independently, emits structured logs (ISSL records) at defined checkpoints, and receives signals from an orchestrator. The orchestrator reads a configuration file that declares which models exist and how they are connected, manages a global simulation clock, routes inter-compartmental signals through optional transfer models, and emits its own unified output log (OISSL) that can be consumed by analysis scripts or a renderer.

This repository provides:

- The **ISSL and OISSL JSON schemas** (in `schemas/`)
- A **reference orchestrator implementation** (in `orchestrator/`)
- **Four model implementations** used in the paper (in `models/`), all calibrated to *Mus musculus* C57BL/6 data
- **Ready-to-run configuration graphs** for all five simulation runs from the paper (in `configs/`)
- **A figure generation script** that renders all paper figures from ISSL logs (in `figures/`)
- **Analysis notebooks** for reproducing paper figures (in `analysis/`)

**Visualisation** (3D Godot renderer) is handled in a separate repository, out of the scope of this repo

---

## Repository structure

```
oisa/
│
├── schemas/                        # JSON schemas and specifications
│   ├── issl_v1.schema.json         # ISSL model log format
│   ├── oissl_v1.schema.json        # OISSL orchestrator output log format
│   ├── config_graph_v1.schema.json # Configuration graph schema
│   └── scene_schemas/
│       └── immune_ontogeny_v1.json # Scene schema for the Godot renderer
│
├── orchestrator/                   # Reference orchestrator implementation
│   ├── README.md                   # Orchestrator-specific documentation
│   ├── main.py                     # Entry point
│   ├── components/
│   │   ├── ingestion.py            # Component 1: ISSL ingestion + validation
│   │   ├── scheduler.py            # Component 2: Global clock + temporal scheduling
│   │   ├── causal_resolver.py      # Component 3: DAG + signal routing
│   │   ├── constraint_engine.py    # Component 4: Biological plausibility checks
│   │   ├── state_registry.py       # Component 5: Global immune state (GIS)
│   │   ├── transfer_dispatcher.py  # Component 6: Transfer model invocation
│   │   ├── calibration_bridge.py   # Component 7: EHR / data ingestion
│   │   ├── output_aggregator.py    # Component 8: OISSL emission
│   │   └── watchdog.py             # Component 9: Health monitoring
│   └── tests/
│
├── models/
│   ├── bm_haematopoiesis/          # Model 1: Bone marrow ODE
│   │   ├── README.md
│   │   ├── model.py                # ODE implementation + ISSL emitter
│   │   ├── parameters.yaml         # Kinetic parameters with sources
│   │   └── tests/
│   │
│   ├── blood_transit/              # Model 2: Blood transit ODE (transfer model)
│   │   ├── README.md
│   │   ├── model.py
│   │   ├── parameters.yaml
│   │   └── tests/
│   │
│   ├── thymus_selection/           # Model 3: Thymic selection ABM
│   │   ├── README.md
│   │   ├── model.py                # ABM implementation + ISSL emitter
│   │   ├── parameters.yaml
│   │   └── tests/
│   │
│   └── peripheral_ln/              # Model 4: Peripheral LN ODE
│       ├── README.md
│       ├── model.py
│       ├── parameters.yaml
│       └── tests/
│
├── configs/                        # Ready-to-run configuration graphs
│   ├── run_BM1_baseline.yaml       # Single bone marrow ODE
│   ├── run_THY1_baseline.yaml      # Single thymus ABM
│   ├── run_COMP1_direct.yaml       # ODE → ABM, no transfer model
│   ├── run_COMP2_transfer.yaml     # ODE → transfer ODE → ABM
│   └── run_COMP3_full_graph.yaml   # All four models
│
├── analysis/
│   ├── notebooks/
│   │   ├── fig2_bm_baseline.ipynb
│   │   ├── fig3_thymus_baseline.ipynb
│   │   ├── fig4_composition.ipynb
│   │   └── fig5_fig6_full_graph.ipynb
│   └── utils/
│       └── oissl_parser.py         # Helper: load + query OISSL logs
│
├── docs/
│   ├── issl_specification.md       # Detailed ISSL field reference
│   ├── config_graph_reference.md   # Configuration graph field reference
│   └── adding_a_model.md           # Step-by-step guide for contributors
│
├── figures/                        # Figure generation
│   ├── generate_figures.py         # Renders all paper figures from ISSL logs
│   ├── BM1_baseline.png
│   ├── THY1_baseline.png
│   ├── COMP1_direct.png
│   ├── COMP2_transfer.png
│   └── COMP3_full_graph.png
│
├── supplementary/
│   └── table_S1_parameters.csv     # All model parameter values + sources
│
├── environment.yml                 # Conda environment
├── pyproject.toml
└── README.md                       # This file
```

---

## Quick start

```bash
# 1. Clone and set up the environment
git clone https://github.com/[org]/oisa.git
cd oisa
conda env create -f environment.yml
conda activate oisa

# 2. Run the full four-model composition (COMP-3 from the paper; 365-day simulation)
PYTHONPATH=. python orchestrator/main.py \
    --config configs/run_COMP3_full_graph.yaml \
    --run-id latest

# 3. Inspect the output logs
ls results/COMP3/latest/          # per-model ISSL logs + orchestrator OISSL

# 4. Generate figures
python figures/generate_figures.py

# 5. Open the analysis notebook for Figure 5
jupyter lab analysis/notebooks/fig5_fig6_full_graph.ipynb
```

---

## Architecture overview

The three components of OISA interact as shown below. Models run as independent processes and communicate with the orchestrator exclusively through ISSL log records. The orchestrator never touches model internals.

```
┌─────────────────────────────────────────────────────────────────┐
│                      Configuration graph                        │
│              (YAML — nodes, edges, clock, renderer)             │
└───────────────────────────┬─────────────────────────────────────┘
                            │ read at init
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│                         Orchestrator                            │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌────────────────┐  │
│  │Scheduler │  │ Causal   │  │Constraint│  │ Output         │  │
│  │(GSimT)   │  │ resolver │  │ engine   │  │ aggregator     │  │
│  └────┬─────┘  └────┬─────┘  └────┬─────┘  └───────┬────────┘  │
│       │              │              │                │           │
│  ┌────┴─────┐  ┌─────┴────┐  ┌─────┴──────┐  ┌─────┴────────┐  │
│  │ISSL      │  │Transfer  │  │State       │  │ Watchdog     │  │
│  │ingestion │  │dispatcher│  │registry    │  │ monitor      │  │
│  └────┬─────┘  └──────────┘  └────────────┘  └──────────────┘  │
└───────┼─────────────────────────────────────────────────────────┘
        │ ISSL records                      OISSL render stream │
        │ (step commands ↓ / logs ↑)                           ▼
   ┌────┴──────┐  ┌───────────┐  ┌──────────┐         ┌───────────┐
   │  Model 1  │  │  Model 2  │  │  Model 3 │         │ Renderer  │
   │  (ODE)    │  │  (ODE)    │  │  (ABM)   │         │ / analysis│
   └───────────┘  └───────────┘  └──────────┘         └───────────┘
```

### The ISSL format

Every model emits an ISSL (Internal Simulation State Log) record at each checkpoint. It is a JSON-LD document with six sections:

| Section | What it contains |
|---|---|
| `envelope` | Model identity, version, `sim_time_s`, schema URI |
| `continuous_state` | Running entity populations: count, unit, fitness, surface markers, `ci_95` |
| `discrete_events` | Punctual biological events: cell death bursts, selection events, cytokine peaks |
| `export_signals` | Inter-compartmental fluxes available for routing to other models |
| `internal_parameters` | Kinetic parameter values with posterior distributions |
| `watchdog` | Health status, OOD flag, divergence score, next checkpoint time |

Entity identifiers use OBO ontology URIs where available (`CL:0000037` for HSC, `GO:0045058` for positive T-cell selection) or namespaced custom identifiers for entities without a current OBO term (e.g. `MODEL_BM_custom:stromal_niche_type_A`).

Full field reference: [`docs/issl_specification.md`](docs/issl_specification.md)
JSON schema: [`schemas/issl_v1.schema.json`](schemas/issl_v1.schema.json)

### The configuration graph

The configuration graph is a YAML file with five top-level keys:

```yaml
oisa_version: "1.0"

models:           # list of model nodes (id, formalism, executable, issl_port, delta_t_s)
                  #   optional: model_args (list of extra CLI arguments passed to the model process)
edges:            # directed signal flows (source, signal_id, target, lag, activation_threshold)
transfer_models:  # optional models on edges (invoked by the transfer dispatcher)
global_clock:     # GSimT parameters (start_s, end_s, checkpoint_interval_s)
renderer:         # OISSL output configuration (format, target, emit_interval_s, scene_schema_uri)
```

The `lag` field on an edge is either `constant:N` (fixed N-second delay) or `model:ID` (the lag is read from the output of a named transfer model). The optional `model_args` list on a model node passes additional command-line arguments to the model subprocess at startup (e.g. `["--baseline-import", "0"]` to disable synthetic thymic import when running in composition). See [`docs/config_graph_reference.md`](docs/config_graph_reference.md) for the full schema.

### The orchestrator

The orchestrator is a nine-component Python server process. It reads the configuration graph at startup, establishes socket connections to all model processes, and manages the simulation run. Its nine components and their responsibilities are documented in [`orchestrator/README.md`](orchestrator/README.md).

The orchestrator emits its own log format — the **OISSL** (Orchestrator ISSL) — at each global checkpoint. The OISSL merges all model state snapshots into a single `global_immune_state` section, records signal routing events in `composition_events`, and includes a `render_entities` section that maps biological entities to 3D scene objects for consumption by the Godot renderer.

---

## Models

### Model 1 — Bone marrow haematopoiesis (ODE)

**Role in the composition:** source model — produces T-lineage progenitor cells exported to the blood transit model or directly to the thymus.

**Formalism:** System of ordinary differential equations. Population-level mean-field dynamics are appropriate here because haematopoietic stem cell (HSC) pool behaviour is well characterised at the population scale and individual stochasticity is not the primary driver of progenitor output variability.

**Biological scope:** Models the bone marrow niche from HSC self-renewal through multipotent progenitor commitment to early T-lineage progenitor (DN1 / ETP) export. Does not model myeloid or B-cell lineages explicitly; they are included as a competing-fate sink term.

**Model version:** v4 (5-compartment). Species: *Mus musculus* C57BL/6, wild-type adult (8–12 weeks).

**Step-by-step model construction:**

1. **Define the state vector.** The model tracks five compartments along the HSC→T-lineage axis: haematopoietic stem cell (HSC), multipotent progenitor (MPP), lymphoid-primed multipotent progenitor (LMPP), common lymphoid progenitor (CLP), and early T-lineage progenitor (DN1/ETP). Each compartment is a continuous scalar variable (cells).

2. **HSC self-renewal and loss.** HSCs self-renew with logistic density dependence toward a niche carrying capacity (`K_niche = 11 000 cells`) and are lost by differentiation into MPPs (`d_HSC_MPP = 0.003 day⁻¹`) and apoptosis (`d_apop = 0.0005 day⁻¹`). At equilibrium, `HSC_eq ≈ 9 075 cells`.

3. **MPP transit amplification.** MPPs undergo transit amplification before lineage commitment: they self-renew at `r_MPP = 0.20 day⁻¹` while exiting to LMPP at `d_MPP_LMPP = 0.15 day⁻¹` and dying at `d_MPP_death = 0.07 day⁻¹`. The net loss rate `λ_MPP = d_MPP_LMPP + d_MPP_death − r_MPP = 0.020 day⁻¹ > 0` ensures stability. At equilibrium, `MPP_eq ≈ 1 361 cells`.

4. **LMPP lymphoid restriction.** A fraction `f_lymphoid = 0.35` of MPP differentiation output is committed to the lymphoid fate (→ LMPP); the remainder enters the myeloid/erythroid sink (not tracked). LMPPs exit to CLP at `d_LMPP_CLP = 0.08 day⁻¹`. At equilibrium, `LMPP_eq ≈ 595 cells`.

5. **CLP and DN1 production.** CLPs are produced from LMPP differentiation and lost to T-lineage commitment (`alpha_T = 0.08 day⁻¹`), B/NK lineages (`d_CLP_other = 0.05 day⁻¹`), and apoptosis. DN1 cells are produced from CLP T-lineage commitment and exit the bone marrow at `export_rate = 0.15 day⁻¹`. At equilibrium, the export flux is approximately **27.2 cells·day⁻¹** — within the published murine range of 10–100 cells·day⁻¹.

6. **Niche signalling feedback (optional).** A feedback term captures the effect of stromal cell-derived cytokines (CXCL12, SCF) on HSC self-renewal rate. In the baseline configuration, this is held at a constant basal level.

7. **ISSL emission.** At each checkpoint (`delta_t_s = 21600`, i.e. every 6 hours), the model serialises all five compartment values with first-order relative-σ `ci_95` (propagated from parameter posterior variances), the surface marker profiles of each population, and the current parameter estimates. The `export_signals` section contains the DN1 export flux with entity ID `CL:0002420` (signal ID: `bm_haematopoiesis.progenitor_export`). These identifiers are stable across model versions and are required for orchestrator signal routing.

8. **Parameter calibration.** All kinetic parameters are drawn from published murine haematopoiesis literature (Busch et al. 2015, Adolfsson et al. 2005, Bhandoola et al. 2007, and others — see `models/bm_haematopoiesis/parameters.yaml`). CI-95 bounds use first-order sensitivity propagation; for MPP in the near-balanced regime (true sensitivity S(r_MPP) ≈ 10), a conservative S=1 bound is reported as a documented lower bound on uncertainty.

**Key outputs logged in ISSL:** HSC, MPP, LMPP, CLP, DN1 counts; DN1 export flux (cells·day⁻¹); `ci_95` on all compartments; parameter posteriors.

---

### Model 2 — Blood transit kinetics (ODE)

**Role in the composition:** transfer model — sits on the edge between Model 1 and Model 3. It is not a biological compartment in its own right but a model of the transit process between two compartments. The orchestrator invokes it when a signal arrives on the `BM.progenitor_export` edge and uses its output to determine (a) how many cells arrive at the thymus and (b) after what delay.

**Formalism:** ODE. Transit kinetics are well described by first-order exponential loss, making ODEs a natural choice.

**Biological scope:** Models the fate of DN1 progenitor cells after egress from the bone marrow sinusoids into the peripheral blood, until their arrival at thymic blood vessels and homing into the thymic cortex.

**Step-by-step model construction:**

1. **Define the transit pool.** The state is a single scalar: the count of DN1 cells currently in transit in the bloodstream.

2. **Input signal.** The model receives the bone marrow `export_signals.flux` value as its initial condition at `t = 0` of the transit simulation. This is the count of cells entering the blood at the moment of BM export.

3. **First-order cell loss.** Cells leave the transit pool by two processes: thymic homing (cells successfully arriving at the thymus) and death/margination (cells that are lost before reaching the thymus). Both are modelled as first-order processes with rate constants calibrated from published data on DN1 progenitor half-life in peripheral blood.

4. **Transit duration.** The characteristic transit time (mean time from BM egress to thymic arrival) is a model parameter. The ODE is integrated until the cumulative homing flux reaches 82% of its asymptotic value (the 82nd percentile, calibrated to the published murine transit window of 1–4 days; Schwarz & Bhandoola 2004, Scimone et al. 2006), yielding a lag of approximately **4 days**. The elapsed integration time is the lag value written to the ISSL `export_signals` section.

5. **ISSL emission.** The model emits a single ISSL record at the end of the transit simulation. The `export_signals` section contains: the number of viable cells delivered to the thymus (`flux`), the computed transit lag in seconds (`lag_s`, which the orchestrator reads to schedule delivery), and the `ci_95` on both quantities.

6. **Parameter calibration.** Transit kinetic parameters are drawn from published data on murine thymocyte precursor trafficking. The death rate accounts for known vulnerability of DN1 cells to oxidative stress during transit.

**Key outputs logged in ISSL:** Cells delivered, transit lag (seconds), cell viability fraction, `ci_95` on all three.

---

### Model 3 — Thymic T-cell selection (ABM)

**Role in the composition:** intermediate model — receives progenitor cells from the blood transit model (or directly from bone marrow in COMP-1) and exports mature naïve T cells to the peripheral lymph node model.

**Formalism:** Agent-based model. Thymic selection is driven by individual stochastic interactions between thymocyte T-cell receptors (TCRs) and self-peptide/MHC complexes on cortical thymic epithelial cells (cTECs) and medullary TECs (mTECs). The repertoire-shaping outcome depends on the distribution of individual binding affinities — a process that cannot be faithfully captured by mean-field ODEs.

**Biological scope:** Models the thymic cortex and medulla as two spatial zones. Tracks the progression of thymocytes through the double-negative (DN) → double-positive (DP) → single-positive (SP) developmental stages, including positive selection in the cortex and negative selection (clonal deletion) in the medulla.

**Step-by-step model construction:**

1. **Define the agent.** Each thymocyte is an agent with the following attributes: developmental stage (DN1, DN2, DN3, DN4, DP, CD4-SP, CD8-SP), a TCR affinity score drawn from a distribution at agent creation, position in the cortex or medulla, age (time steps since entry), and a fate status (alive, positively selected, negatively selected / deleted, exported).

2. **Spatial environment.** The thymus is represented as two zones: a cortical zone and a medullary zone. In the cortex, each DP thymocyte encounters cTEC-presented self-peptide/MHC at each time step; in the medulla, SP thymocytes encounter mTEC/dendritic cell-presented tissue-restricted antigens. Encounter probability per time step is a parameter.

3. **TCR affinity distribution.** At agent creation, each thymocyte is assigned a TCR affinity score drawn from a log-normal distribution. This distribution encodes the diversity of the naïve T-cell repertoire entering the thymus.

4. **Positive selection rule.** At each cortical encounter, a DP thymocyte undergoes positive selection if its TCR affinity score falls within the window [θ_low, θ_high]. Cells with affinity below θ_low die by neglect; cells above θ_high proceed to negative selection assessment. The two thresholds are model parameters.

5. **Negative selection rule.** In the medulla, positively selected SP thymocytes encounter tissue-restricted self antigens. If the TCR affinity for any encountered antigen exceeds θ_negative, the cell undergoes clonal deletion (apoptosis). The probability of encountering any given antigen per medullary dwell time follows a Poisson distribution parameterised by the medullary antigen diversity and encounter rate.

6. **Migration rules.** Agents migrate between cortex and medulla at age-dependent rates reflecting the known temporal progression of thymic education. DN → DP transition occurs in the cortex; SP cells migrate to the medulla after positive selection. Export to the periphery occurs after medullary dwell time completion without deletion.

7. **Import signal handling.** At each checkpoint, the model queries its Accept function for pending import signals from the orchestrator. New DN1 agent instances are created with TCR affinities drawn from the distribution and placed in the cortical DN1 zone. The number of new agents is sampled from the `export_signals.flux` ± `ci_95` of the incoming signal.

8. **ISSL emission.** At each 24-hour checkpoint, the model aggregates agent counts by developmental stage and zone into the `continuous_state` section. Individual selection events (positive selection completions, clonal deletions) are recorded in the `discrete_events` section with counts, `n_realisations`, and variance across ensemble runs. The `export_signals` section contains the daily naïve T cell export flux with entity ID `CL:0000898`.

9. **Stochastic variance.** The model is run for `n_realisations = 12` independent realisations per checkpoint (configurable). The `ci_95` fields in the ISSL are computed empirically across realisations. The orchestrator reconciles this stochastic output with any downstream ODE model by using the distributional mean as the point input.

10. **Parameter calibration.** TCR affinity thresholds, encounter probabilities, and migration rates are calibrated against published data on positive/negative selection yields in wild-type murine thymus (see `supplementary/table_S1_parameters.csv`).

**Key outputs logged in ISSL:** Agent counts per developmental stage, daily selection event counts (positive / negative / neglect death), naïve T export flux, ensemble variance, parameter posteriors.

---

### Model 4 — Peripheral lymph node dynamics (ODE)

**Role in the composition:** sink model in the baseline configuration — receives naïve T cells from the thymus and models their subsequent fate in the peripheral lymph node pool. Can also be a source model in extended configurations that include downstream effector or tissue models.

**Formalism:** ODE. The peripheral naïve T-cell pool is maintained by a balance of thymic export, homeostatic proliferation, and peripheral turnover — processes that are well characterised at the population level and where mean-field ODEs are appropriate.

**Biological scope:** Models the naïve CD4+ and CD8+ T-cell pools in the peripheral lymph node. Tracks homeostatic proliferation driven by IL-7 signalling, peripheral T-cell loss by activation (transition to effector/memory; modelled as a sink), and natural death. Does not model antigen-specific activation in the current implementation; this is reserved for a future extension.

**Step-by-step model construction:**

1. **Define the state vector.** Two continuous scalar variables: naïve CD4+ T-cell pool size (cells) and naïve CD8+ T-cell pool size (cells). The CD4/CD8 ratio at thymic export (approximately 2:1) is used to partition the incoming signal between the two pools.

2. **Thymic import.** At each time step, the model's Accept function receives the `THY.naive_T_export` signal from the orchestrator (delivered with a fixed 2-day homing lag declared in the configuration graph). The import is added as a positive term to both pool ODEs, partitioned by the CD4/CD8 split ratio.

3. **Homeostatic proliferation.** Each naïve T-cell pool proliferates at a rate proportional to the difference between the current pool size and the homeostatic set point, driven by IL-7 availability. The set point and the proliferation rate constant are model parameters. This term implements the compensatory proliferation observed after lymphopenia.

4. **Peripheral loss.** Naïve T cells are lost by spontaneous activation (a small constant per-cell rate representing tonic TCR signalling and bystander activation) and by natural death. Both are first-order processes.

5. **IL-7 signalling (optional feedback).** An optional input signal from a stromal cell model can modulate the IL-7 level and therefore the homeostatic proliferation rate. In the current COMP-3 configuration, IL-7 is held at a constant basal level.

6. **ISSL emission.** At each 12-hour checkpoint, the model emits pool sizes for CD4+ and CD8+ naïve T cells with analytical `ci_95`, the current import rate (from the last received signal), and the homeostatic proliferation rate. The `export_signals` section is empty in the COMP-3 configuration but would contain effector T-cell flux in an extended model.

7. **Parameter calibration.** Homeostatic parameters are drawn from published data on murine peripheral T-cell pool kinetics following thymic reconstitution.

**Key outputs logged in ISSL:** CD4+ naïve T count, CD8+ naïve T count, CD4/CD8 ratio, homeostatic proliferation rate, import rate, `ci_95` on all variables.

---

## Running the simulations

### Simulation runs from the paper

All five runs from the paper are provided as ready-to-run configuration files:

| Config file | Description | Paper section |
|---|---|---|
| `run_BM1_baseline.yaml` | Bone marrow ODE in isolation, 30 days | §5.1 |
| `run_THY1_baseline.yaml` | Thymus ABM in isolation, constant synthetic input, 30 days | §5.1 |
| `run_COMP1_direct.yaml` | BM-ODE → Thymus-ABM, direct coupling, no lag | §5.2 |
| `run_COMP2_transfer.yaml` | BM-ODE → transit-ODE → Thymus-ABM, model-derived lag (~4 days) | §5.3 |
| `run_COMP3_full_graph.yaml` | All four models, full immune ontogeny graph, **365 days** | §5.4 |

Run any of them with:

```bash
python orchestrator/main.py --config configs/run_COMP3_full_graph.yaml \
    --output-dir logs/COMP3 \
    --log-level INFO
```

Output is written to `logs/COMP3/`:
- `logs/COMP3/issl/BM_haematopoiesis_v4/` — per-checkpoint ISSL records from Model 1
- `logs/COMP3/issl/Thymus_selection_v2/` — per-checkpoint ISSL records from Model 3
- (... etc. for all models)
- `logs/COMP3/oissl/` — per-checkpoint OISSL records from the orchestrator

### Writing your own configuration graph

See [`docs/config_graph_reference.md`](docs/config_graph_reference.md) for the full field reference. The minimal valid configuration for a single model is:

```yaml
oisa_version: "1.0"

models:
  - id: my_model
    formalism: ODE          # or ABM, or hybrid
    executable: "./models/my_model/model.py"
    issl_port: "tcp://localhost:5010"
    delta_t_s: 86400

global_clock:
  start_s: 0
  end_s: 2592000
  checkpoint_interval_s: 86400
```

---

## Output logs and analysis

### ISSL per-model logs

Each model writes one JSON file per checkpoint to `logs/<run_id>/issl/<model_id>/checkpoint_<sim_time_s>.json`. The file conforms to `schemas/issl_v1.schema.json`.

To load and inspect all checkpoints for a single model:

```python
from analysis.utils.oissl_parser import load_issl_series

records = load_issl_series("logs/COMP3/issl/BM_haematopoiesis_v4/")

# records is a list of dicts, one per checkpoint, sorted by sim_time_s
hsc_counts = [r["continuous_state"][0]["count"] for r in records]
```

### OISSL orchestrator logs

The orchestrator writes one OISSL file per global checkpoint to `logs/<run_id>/oissl/checkpoint_<sim_time_s>.json`. Each OISSL record contains:

- `global_immune_state` — merged snapshot of all model states at this checkpoint
- `composition_events` — signal routing events, constraint violations, watchdog alerts
- `render_entities` — entity states mapped to scene objects for the renderer

```python
from analysis.utils.oissl_parser import load_oissl_series

oissl = load_oissl_series("logs/COMP3/oissl/")

# End-to-end latency: time from BM export signal to peripheral pool update
bm_export_times  = [e["sim_time_s"] for r in oissl
                    for e in r["composition_events"]
                    if e["type"] == "signal_emitted" and e["signal_id"] == "BM.progenitor_export"]
periph_import_times = [e["sim_time_s"] for r in oissl
                       for e in r["composition_events"]
                       if e["type"] == "signal_delivered" and e["target_model"] == "PeripheralLN_ODE"]
```

### Render stream

If a `renderer` block is present in the configuration graph, the orchestrator emits OISSL records to the specified target (a file path or a TCP socket). The Godot renderer (separate repository) connects to this socket and updates the scene at each checkpoint.

To write the render stream to a file for offline replay:

```yaml
renderer:
  format: "OISA-render-v1"
  target: "file:logs/COMP3/render_stream.ndjson"
  emit_interval_s: 86400
  scene_schema_uri: "schemas/scene_schemas/immune_ontogeny_v1.json"
```

---

## Extending OISA

### Adding a new model

1. Create a directory under `models/your_model_name/`.
2. Implement `model.py` with two required functions:
   - `step(sim_time_s, incoming_signals)` — advance the model by one `delta_t_s` given any pending input signals; return nothing.
   - `emit_issl(sim_time_s)` — serialise the current state as a dict conforming to `schemas/issl_v1.schema.json`; return the dict.
3. Write `parameters.yaml` listing all kinetic parameters with values, units, and literature sources.
4. Register the model in a configuration graph YAML.

Detailed walkthrough: [`docs/adding_a_model.md`](docs/adding_a_model.md).

### Adding a transfer model on an edge

A transfer model is any OISA-compliant model whose `export_signals` section contains a field named `lag_s` (the computed delay in seconds). The orchestrator reads this field automatically when the edge `lag` is set to `model:ID`.

```yaml
edges:
  - source: my_source_model
    signal_id: my_source_model.export_flux
    target: my_target_model
    lag: "model:my_transfer_model"     # orchestrator invokes this model on the edge
```

The transfer model receives the source export signal as its initial condition (via `Accept`), runs its own simulation, and emits an ISSL record with the transformed signal in `export_signals` and the computed lag in `export_signals[0].lag_s`.

### Using custom (non-ontologised) entities

If your model tracks entities that do not have an OBO term, use `entity_class: "custom"` in the ISSL record:

```json
{
  "entity_class": "custom",
  "entity_id": "MODEL_BM_custom:stromal_niche_type_A",
  "label": "CXCL12-abundant reticular cell",
  "count": 1240,
  "unit": "cells"
}
```

The orchestrator will track this entity in the state registry and pass it through the render stream. It will not apply ontology-based constraint checks to custom entities. Custom entity namespaces should be declared in the `wildcard_namespace` block of the configuration graph to enable future mapping to OBO terms.

---

## Schemas and specifications

| File | Description |
|---|---|
| `schemas/issl_v1.schema.json` | JSON Schema for ISSL model log records |
| `schemas/oissl_v1.schema.json` | JSON Schema for OISSL orchestrator output records |
| `schemas/config_graph_v1.schema.json` | JSON Schema for configuration graph files |
| `schemas/scene_schemas/immune_ontogeny_v1.json` | Scene schema for the Godot immune ontogeny renderer |
| `docs/issl_specification.md` | Full narrative field reference for ISSL |
| `docs/config_graph_reference.md` | Full narrative field reference for the configuration graph |

---

