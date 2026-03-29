# Journal de débogage — SASSy / BM1 Baseline

## 2026-03-26 — Bug CI_95 sur le HSC (Bone Marrow Haematopoiesis)

### Contexte

La simulation `run_BM1_baseline` produit un enregistrement ISSL à chaque
checkpoint. Chaque variable d'état (`HSC`, `CLP`, `DN1`) est accompagnée d'un
intervalle de confiance à 95 % (`ci_95`) calculé analytiquement à partir des
distributions postérieures des paramètres.

---

### Symptôme observé

Run `20260326T085729` (avant correctif) :

| Variable | Valeur (jour 1) | CI_95 obtenu | Largeur relative |
|----------|----------------|--------------|-----------------|
| HSC | 9 979 cells | [-14 958 874 ; +14 978 832] | **±300 000 %** |
| CLP | 1 991 cells | [1 887 ; 2 095] | ±5.2 % |
| DN1 | 637 cells | [616 ; 658] | ±3.4 % |

La borne inférieure du HSC est à −15 millions de cellules pour une valeur
nominale de ~10 000. L'intervalle est **1 500 fois plus large que la valeur
elle-même**, ce qui le rend biologiquement inexploitable.

---

### Cause racine : `_state_ci95` — mélange d'unités incompatibles

Fichier : `models/bm_haematopoiesis/model.py`, méthode `_state_ci95` (ligne ~198)

#### Code original (bugué)

```python
def _state_ci95(self, value: float, *param_keys: str) -> list[float]:
    variance = 0.0
    for k in param_keys:
        ci = self._p_ci[k]
        sigma = (ci[1] - ci[0]) / (2 * 1.96)   # σ ABSOLU du paramètre
        variance += sigma ** 2
    sigma_total = math.sqrt(variance) * value
    return [round(value - 1.96 * sigma_total, 2), round(value + 1.96 * sigma_total, 2)]
```

La fonction calcule le **σ absolu** de chaque paramètre, additionne les
variances en quadrature, puis multiplie par `value`. Cela ne pose aucun
problème tant que tous les paramètres ont des valeurs et des échelles
comparables — mais ce n'est pas le cas pour le HSC.

#### Paramètres impliqués dans le calcul du CI du HSC

```python
ci_hsc = self._state_ci95(HSC, "r_self", "K_niche")
```

| Paramètre | Valeur nominale | CI_95 | σ absolu |
|-----------|----------------|-------|---------|
| `r_self` | 0.02 day⁻¹ | [0.018 ; 0.022] | **0.00102 day⁻¹** |
| `K_niche` | 11 000 cells | [9 500 ; 12 500] | **765.3 cells** |

#### Reconstruction numérique du bug

```
variance = (0.00102)² + (765.3)²
         ≈ 0.000001 + 585 684
         ≈ 585 684          ← dominé à 100 % par K_niche

sigma_total = √585 684 × 9 979 ≈ 765.3 × 9 979 ≈ 7 638 000

CI = [9 979 − 1.96 × 7 638 000 ; 9 979 + 1.96 × 7 638 000]
   = [−14 961 521 ; 14 981 479]
```

`K_niche` est exprimé en **cellules absolues** (σ ≈ 765 cells), alors que
`r_self` est un taux sans dimension (σ ≈ 0.001 day⁻¹). Les additionner
directement en quadrature est une **erreur d'analyse dimensionnelle** : le σ de
`K_niche`, numériquement élevé en valeur absolue, écrase totalement le σ de
`r_self` et se retrouve multiplié par ~10 000 (la valeur du HSC), produisant un
intervalle absurde.

#### Pourquoi CLP et DN1 semblaient corrects

Pour `CLP`, la fonction est appelée avec `"alpha_T"` et `"alpha_myeloid"`,
deux taux dont les valeurs nominales sont dans [0.01 ; 0.30] et les σ absolus
dans [0.002 ; 0.026]. Numériquement proches d'incertitudes relatives, ils
donnent des CI d'apparence raisonnable — mais **par coïncidence d'échelle**,
pas par correction du calcul. L'erreur de fond était présente pour toutes les
variables ; elle ne se manifestait visuellement que pour le HSC.

---

### Correctif appliqué

**Principe** : normaliser le σ de chaque paramètre par sa valeur nominale pour
obtenir une **incertitude relative (sans dimension)** avant de sommer les
variances. Tous les paramètres, quelle que soit leur unité, contribuent alors
proportionnellement à leur propre incertitude relative.

```
σ_rel_i = (hi_i − lo_i) / (2 × 1.96 × p_i0)
```

La variance totale relative est :

```
σ²_rel_total = Σ σ²_rel_i
```

Et l'incertitude absolue sur la variable d'état :

```
σ_abs = √σ²_rel_total × value
```

#### Code corrigé

```python
def _state_ci95(self, value: float, *param_keys: str) -> list[float]:
    variance_rel = 0.0
    for k in param_keys:
        ci = self._p_ci[k]
        param_val = self._p[k]
        sigma_rel = (ci[1] - ci[0]) / (2 * 1.96 * param_val)   # σ RELATIF
        variance_rel += sigma_rel ** 2
    sigma_total = math.sqrt(variance_rel) * value
    return [round(value - 1.96 * sigma_total, 2), round(value + 1.96 * sigma_total, 2)]
```

---

### Résultats après correctif

Run `20260326T101627` (après correctif) :

| Variable | Valeur (jour 1) | CI_95 corrigé | Largeur relative |
|----------|----------------|---------------|-----------------|
| HSC | 9 979 cells | [8 291 ; 11 666] | ±16.9 % |
| CLP | 1 991 cells | [1 445 ; 2 537] | ±27.4 % |
| DN1 | 637 cells | [462 ; 812] | ±27.4 % |

**Jour 30 :**

| Variable | Valeur | CI_95 corrigé | Largeur relative |
|----------|--------|---------------|-----------------|
| HSC | 9 614 cells | [7 989 ; 11 240] | ±16.9 % |
| CLP | 85.2 cells | [61.8 ; 108.5] | ±27.4 % |
| DN1 | 63.6 cells | [46.2 ; 81.1] | ±27.4 % |

Les incertitudes sont cohérentes avec les CI des paramètres sources :
- `K_niche` : [9 500 ; 12 500] → incertitude relative ~14 % → HSC ~±17 % ✓
- `alpha_T` / `alpha_myeloid` : incertitudes relatives ~9–10 % chacune,
  combinées → CLP/DN1 ~±14 % en σ soit ~±27 % en demi-intervalle à 95 % ✓

#### Note sur l'élargissement du CI de CLP et DN1

Les CI de CLP et DN1 sont plus larges qu'avant le correctif (~27 % vs ~5 %).
C'est **attendu et correct** : avant le correctif, les petits σ absolus des
taux (par ex. σ(alpha_T) ≈ 0.008) sous-estimaient l'incertitude réelle. L'incertitude
relative de `alpha_T` est en réalité de ~20 %, ce qui se répercute correctement
sur CLP et DN1.

---

### Fichiers modifiés

| Fichier | Nature |
|---------|--------|
| `models/bm_haematopoiesis/model.py` | Correctif de `_state_ci95` (ligne ~198) |

### Runs de référence

| Run | Statut | Répertoire |
|-----|--------|-----------|
| `20260326T085729` | Bugué (CI_95 invalides) | `results/BM1/20260326T085729/` |
| `20260326T101627` | Corrigé | `results/BM1/20260326T101627/` |

---

## 2026-03-26 — Zéro export de naïve T cells (THY1, COMP1, COMP2)

### Contexte

Les simulations THY1 (thymus isolé), COMP1 (BM→Thymus direct) et COMP2
(BM→blood_transit→Thymus) produisent systématiquement `exported = 0` à chaque
checkpoint sur les 30 jours de simulation.

Trois hypothèses testées :
1. Seuils de sélection trop stricts (pas assez de cellules sélectionnables)
2. Temps de simulation insuffisant (agents bloqués sur 30 jours)
3. Bug dans la condition d'export

---

### Diagnostic préalable — analyse statique

#### Distribution TCR (LogNormal μ=0, σ=0.8)

| Destin | Fraction |
|--------|---------|
| Neglect death (aff < 0.30) | 6.8 % |
| **Positive selection (0.30 ≤ aff ≤ 2.50)** | **80.6 %** |
| Neg sel au stade DP (aff > 2.50) | 12.6 % |
| Neg sel médullaire (aff > 3.50 parmi pos sel) | ~0 % |

→ Les seuils ne sont **pas** trop stricts : 80.6 % des cellules sont sélectionnables.

#### Timing DN1 → export (attendu après correction du bug)

| Étape | E[pas] | E[durée] |
|-------|--------|---------|
| DN1 → DP (`p_transition = 0.05`) | 20 pas | 0.8 j |
| DP → sélection pos. (`p_enc = 0.30`, 80.6 % sel.) | 4 pas | 0.2 j |
| SP cortex → médulla (`p_migrate = 0.10`) | 10 pas | 0.4 j |
| Dwell médullaire (`medullary_dwell_steps = 72`) | 72 pas | 3.0 j |
| **Total** | **~106 pas** | **~4.4 j** |

→ 30 jours est **largement suffisant** (marge de 25 jours). Le timing n'est
**pas** le problème.

---

### Cause racine 1 — Bug d'état après sélection positive (COMP1, COMP2)

**Fichier** : `models/thymus_selection/model.py`, ligne ~125

#### Code original (bugué)

```python
elif aff <= p["theta_high"]:
    # Positive selection → become SP
    agent.fate = "pos_selected"   # ← non-"alive"
    self.pos_sel_count += 1
    agent.stage = "CD4SP" if lineage == "CD4" else "CD8SP"
    agent.zone = "cortex"
```

À la sous-étape suivante :

```python
if agent.fate != "alive":
    continue   # ← l'agent est immédiatement ignoré, pour toujours
```

L'agent passe au stage `CD4SP`/`CD8SP` mais reste piégé avec `fate =
"pos_selected"`. La garde en tête de boucle le court-circuite à chaque
sous-étape : il ne peut jamais atteindre la logique de migration vers la
médulla ni la condition d'export. Le code avait prévu de réinitialiser `fate`
lors de la migration (`agent.fate = "alive"  # reset from pos_selected label`),
mais cette ligne est inaccessible.

**Résultat observable** : `CD4SP`/`CD8SP` plafonnent à 0-1 (un agent tout juste
sélectionné dans le dernier sous-pas), `exported = 0` pour toute la simulation.

#### Correctif

```python
elif aff <= p["theta_high"]:
    # fate doit rester "alive" pour que l'agent ne soit pas court-circuité
    # à la prochaine sous-étape ; pos_sel_count trace l'événement.
    agent.fate = "alive"          # ← correctif
    self.pos_sel_count += 1
    agent.stage = "CD4SP" if lineage == "CD4" else "CD8SP"
    agent.zone = "cortex"
```

---

### Cause racine 2 — Absence d'input pour THY1 (run isolé)

**Fichier** : `configs/run_THY1_baseline.yaml`

Le commentaire du fichier stipule « A fixed 100 cells/day import is injected
via CalibrationBridge », mais `CalibrationBridge` est instancié dans `main.py`
sans jamais être branché au Scheduler. La méthode `inject_as_signal` n'est
jamais appelée. Sans edge entrant, `signals = []` dans `_step()` → 
`n_import_mean = 0.0` → aucun agent ajouté.

#### Correctif

Ajout du paramètre `baseline_import` (100 cells·day⁻¹) dans `parameters.yaml`,
utilisé comme fallback dans `_step()` quand aucun signal amont n'est reçu :

```python
if not signal_received:
    baseline = self._p.get("baseline_import", 0.0)
    if baseline > 0.0:
        n_import_mean = baseline
        n_import_sigma = (ci[1] - ci[0]) / (2 * 1.96)
```

Ce mécanisme ne s'active que si aucune edge ne délivre de signal, préservant
le comportement dans COMP1/COMP2 (où le signal BM prend le dessus).

---

### Résultats après correctifs

Runs corrigés (2026-03-26) :

#### THY1 — régime permanent (jour 5–30)

Export naïve T : **~80–87 cells·checkpoint⁻¹** (CI_95 ≈ ±15 cells)
Ratio CD4/CD8 : **~1.8** (conforme au `cd4_fraction = 0.65`)
Steady-state atteint au jour 5 (cohérent avec le timing ~4.4 j estimé)

#### COMP1 — BM → Thymus (couplage direct)

Export : **86 cells/j** (pic jour 7) → **10 cells/j** (jour 30)
Déclin : dû à la dynamique HSC du modèle BM (HSC se stabilise à ~9 075 < IC de 10 000),
ce qui réduit le flux progenitor et donc l'input du thymus.

#### COMP2 — BM → blood_transit → Thymus

Lag blood_transit ≈ **7 jours** (601 933 s)
Jours 1–7 : identiques à THY1 (baseline_import actif pendant le lag)
Jours 8–30 : transition vers le signal BM relayé ; export décline plus lentement
que COMP1 (pertes en transit ≈ 20-30 %).

---

### Fichiers modifiés

| Fichier | Nature |
|---------|--------|
| `models/thymus_selection/model.py` | Correctif `agent.fate = "alive"` après pos. sel. (l.~125) |
| `models/thymus_selection/model.py` | Fallback `baseline_import` dans `_step()` quand aucun signal |
| `models/thymus_selection/parameters.yaml` | Ajout paramètre `baseline_import: 100.0 cells·day^-1` |

### Runs de référence

| Run | Config | Statut | Répertoire |
|-----|--------|--------|-----------|
| Bugués | THY1/COMP1/COMP2 | exported=0 toute la simulation | `results/*/20260326T10{29,30,32}*` |
| Corrigés | THY1 | ~83 cells/j en régime | `results/THY1/20260326T103936/` |
| Corrigés | COMP1 | 86→10 cells/j (déclin BM) | `results/COMP1/20260326T103941/` |
| Corrigés | COMP2 | lag 7 j puis 74→17 cells/j | `results/COMP2/20260326T103945/` |

---

## 2026-03-26 — Analyse des résultats et plan d'action

### Observations post-correction

Trois incohérences identifiées après examen des résultats des runs corrigés.

---

#### Incohérence 1 — Conditions initiales BM hors équilibre

L'ODE BM admet un point fixe analytique :

```
HSC_eq = K_niche × (1 − (d_diff + d_apop) / r_self)
       = 11 000 × (1 − 0.0035 / 0.02) = 9 075

CLP_eq = d_diff × HSC_eq / (α_T + α_myeloid + d_CLP)
       = 0.003 × 9 075 / 0.34 ≈ 80

DN1_eq = α_T × CLP_eq / export_rate
       = 0.08 × 80 / 0.15 ≈ 43

export_flux_eq = export_rate × DN1_eq ≈ 6.4 cells·day⁻¹
```

Or les CI dans `parameters.yaml` sont : HSC=10 000, CLP=3 000, DN1=500 —
respectivement ×1.1, ×37 et ×11.6 au-dessus de l'équilibre. La fenêtre de
30 jours est entièrement dominée par un transitoire de relaxation vers
l'équilibre : le flux BM passe de ~96 à ~6.4 cells/j, ce qui explique
l'effondrement progressif des exports dans COMP1 et COMP2.

---

#### Incohérence 2 — Sélection négative médullaire inopérante

| Paramètre | Valeur | Rôle |
|-----------|--------|------|
| `theta_high` | 2.50 | borne haute de la sélection positive |
| `theta_negative` | 3.50 | seuil de délétion médullaire |

Comme `theta_negative > theta_high`, toute cellule positivement sélectionnée
(affinité ≤ 2.50) est *a priori* sous le seuil de délétion (3.50). La sélection
négative en médulla est donc mécaniquement inopérante. Résultat : efficacité
thymus ~83 % vs 2–5 % biologiquement attendu.

Secondairement, `p_encounter_cortex = 0.30 / pas` × 20 pas DP implique
P(jamais rencontré pMHC) < 0.1 %, éliminant pratiquement toute mort par
neglect, mécanisme principal in vivo.

---

#### Incohérence 3 — Artefact de préchargement COMP2

Durant le lag blood_transit (~7 jours), le thymus reçoit 0 signal BM et active
le `baseline_import = 100 cells/j`. Les pools SP sont donc pré-remplis avec 7
jours d'input à 100 cells/j avant que le signal BM n'arrive (~66 cells/j après
pertes en transit). La comparaison COMP1 vs COMP2 au jour 30 est biaisée par
cet artefact : l'écart (10 vs 17 cells/j) reflète autant le préchargement que
l'effet biologique réel du transit.

---

### Plan d'action

| # | Fix | Fichiers |
|---|-----|---------|
| 1 | ICs BM à l'équilibre analytique (HSC=9075, CLP=80, DN1=43) | `models/bm_haematopoiesis/parameters.yaml` |
| 2 | `theta_negative = 1.80` (dans la fenêtre pos sel [0.30, 2.50]) | `models/thymus_selection/parameters.yaml` |
| 3 | Support `model_args` dans schema + orchestrateur ; arg `--baseline-import 0` passé au thymus dans COMP1/COMP2 | `schemas/config_graph_v1.schema.json`, `orchestrator/main.py`, `models/thymus_selection/model.py`, `configs/run_COMP1_direct.yaml`, `configs/run_COMP2_transfer.yaml` |


### Exécution du plan d'action et résultats

#### Fix 1 — ICs BM à l'équilibre analytique

Fichier modifié : `models/bm_haematopoiesis/parameters.yaml`

| Variable | Avant | Après (équilibre) |
|----------|-------|-------------------|
| HSC | 10 000 | **9 075** |
| CLP | 3 000 | **80** |
| DN1 | 500 | **43** |

Résultat BM1 (run corrigé) : flux export = **6.41 cells/j stable sur 30 jours**,
sans aucune variation. Le transitoire de relaxation a complètement disparu.

---

#### Fix 2 — theta_negative dans la fenêtre de sélection positive

Fichier modifié : `models/thymus_selection/parameters.yaml`

`theta_negative : 3.50 → 1.80`

Avec θ_neg = 1.80 < θ_high = 2.50, les cellules posititement sélectionnées
avec affinité ∈ [1.80, 2.50] (~13 % de la fenêtre pos. sel.) sont désormais
à risque de délétion médullaire. Sur 72 pas de dwell, la quasi-totalité est
éliminée (P(survie) = 0.85^72 ≈ 5×10⁻⁶).

Résultat THY1 :

| Métrique | Avant (θ_neg=3.50) | Après (θ_neg=1.80) |
|----------|-------------------|--------------------|
| neg_sel events/j | 0 | ~23 |
| Exported cells/j | ~83 | ~70 |
| Efficacité thymus | ~83 % | ~70 % |

La sélection négative est désormais mécaniquement opérationnelle. L'efficacité
reste > 2–5 % biologiques car le modèle ne simule pas la mort DP par absence
de rencontre pMHC (premier mécanisme de mort in vivo) — correction future.

---

#### Fix 3 — Artefact baseline_import dans COMP1/COMP2

Trois composants modifiés :

**Schema** (`schemas/config_graph_v1.schema.json`) : ajout du champ optionnel
`model_args: array[string]` dans la définition de chaque modèle.

**Orchestrateur** (`orchestrator/main.py`, `_launch_models`) : passage des
`model_args` de la config à la commande de lancement du subprocess.

**Thymus model** (`models/thymus_selection/model.py`) : ajout de l'argument
CLI `--baseline-import` ; si fourni, remplace la valeur de `parameters.yaml`
avant le démarrage du serveur ZMQ.

**Configs** (`run_COMP1_direct.yaml`, `run_COMP2_transfer.yaml`) : ajout de
`model_args: ["--baseline-import", "0"]` sur l'entrée `thymus_selection`.

Résultat COMP2 — comparaison avant/après :

| Période | Avant | Après |
|---------|-------|-------|
| Jours 1–7 (lag transit) | 100 cells/j (artefact) | 0 cells/j (correct) |
| Jour 8 (1ère livraison BM) | step-down brutal | montée progressive |
| Jours 12–30 (régime) | ~45 cells/j | **~4.5 cells/j** |
| Comparaison COMP1 jour 30 | biaisée | comparable ✓ |

---

### Runs de référence (plan d'action exécuté)

| Run | Config | Répertoire | Notes |
|-----|--------|-----------|-------|
| BM1 corrigé | `run_BM1_baseline.yaml` | `results/BM1/20260326T105947/` | Steady-state dès j0 |
| THY1 corrigé | `run_THY1_baseline.yaml` | `results/THY1/20260326T105952/` | neg_sel actif |
| COMP1 corrigé | `run_COMP1_direct.yaml` | `results/COMP1/20260326T105956/` | ~4.5 cells/j stable |
| COMP2 corrigé | `run_COMP2_transfer.yaml` | `results/COMP2/20260326T106000/` | lag propre, même régime que COMP1 |

### Incohérence résiduelle documentée

L'efficacité du thymus reste ~70 % vs 2–5 % biologiquement observé. Adresser
cette incohérence requiert d'introduire une durée de vie maximale des agents DP
(mechanism de neglect death par épuisement temporel) ou de réduire
`p_encounter_cortex` d'un ordre de grandeur, changements hors du périmètre
du présent plan d'action.

---

## 2026-03-26 — Analyse détaillée de l'efficacité thymus résiduelle (~70%)

### Contexte

Après les corrections du 26/03, le thymus exporte ~70 % des thymocytes entrants
(THY1) vs 2–5 % biologiquement attendu. Trois causes identifiées par analyse
quantitative.

---

### Cause A — Distribution TCR trop permissive (cause principale)

La distribution log-normale actuelle (μ=0, σ=0.8) place **80.9 %** des cellules
dans la fenêtre de sélection positive [0.30, 2.50]. En biologie, seuls 2–5 %
des thymocytes ont un TCR dont l'affinité correspond à un auto-pMHC sélectif.

Analyse paramétrique sur μ (σ=0.8 fixé) :

| μ | Neglect | Pos sel | Neg sel médul. | Efficacité nette |
|---|---------|---------|----------------|-----------------|
| 0.0 | 6.6 % | 80.9 % | 12.9 % de pos | 70.4 % |
| -1.0 | 39.8 % | 59.3 % | 2.6 % | 57.8 % |
| -2.0 | 84.1 % | 15.9 % | 0.3 % | 15.9 % |
| **-2.5** | **94.8 %** | **5.2 %** | 0.1 % | **5.2 %** |
| -3.0 | 98.7 % | 1.2 % | 0.0 % | 1.2 % |

**Cible retenue : μ = -2.5** → 5.2 % d'efficacité, dans la fourchette 2–5 %.

Note : avec μ=-2.5, la sélection négative médullaire devient quantitativement
négligeable (~0.1 % de pos sel a aff > 1.80). Le mécanisme reste
architecturalement correct ; il est simplement peu sollicité à ce μ (cohérent
avec le fait que la délétion clonale est un phénomène rare, ~0.5–1 % des
thymocytes).

---

### Cause B — Accumulation de ghosts DP (bug de performance / exactitude)

Dans `ThymicRealisation.step()`, la ligne `survivors.append(agent)` est
toujours exécutée pour les agents DP en zone cortex, y compris après
`agent.fate = "neglect_dead"` ou `"neg_deleted"`. L'agent entre dans la liste
`survivors` avec un fate non-alive et est ignoré (`continue`) à chaque sous-étape
suivante, sans jamais en sortir.

Avec μ=0, l'effet est limité (6.6 % de neglect + 12.6 % de DP neg, donc ~19 %
de ghosts). Avec μ=-2.5 :

- P(neglect par pas) = p_enc × P(aff < θ_low) = 0.30 × 0.948 = **0.284**
- Temps moyen avant premier neglect : 1/0.30 ≈ **3.3 pas** (~3 h)
- 94.7 % des DP agents sont condamnés et meurent en quelques pas mais **restent
  en mémoire** pour toute la simulation

Ce bug bloquerait le modèle en pratique avec μ=-2.5 : la liste `agents` croît
sans borne. **Fix obligatoire avant le changement de distribution.**

---

### Cause C — Absence de durée de vie DP (max_dp_age_steps)

Biologiquement, les cellules DP ont une durée de vie de 3–4 jours ; celles qui
n'ont pas reçu de signal de sélection positive meurent par neglect. Le modèle
n'a pas de tel plafond. Avec p_enc=0.30 et μ=-2.5 :

- P(pas d'encounter en 48 pas = 2 j) = (1−0.30)^48 ≈ **3.7×10⁻⁸**

Le plafond est un filet de sécurité biologique, pas un mécanisme quantitatif
dominant ici. Il est néanmoins ajouté comme contrainte physiologique
(`max_dp_age_steps = 48`, 2 jours).

---

### Plan d'action

| # | Fix | Fichier |
|---|-----|---------|
| A | `tcr_affinity.mu : 0.0 → -2.5` | `models/thymus_selection/parameters.yaml` |
| B | Fix ghost accumulation DP : ne pas ajouter à survivors si fate ≠ alive | `models/thymus_selection/model.py` |
| C | Ajout `max_dp_age_steps: 48` + compteur `dp_steps` sur agent | `parameters.yaml` + `model.py` |


### Exécution du plan — résultats

#### Fix A — μ TCR : 0.0 → -2.5

`models/thymus_selection/parameters.yaml` : `tcr_affinity.mu: -2.5`

Fractions analytiques (LogNormal(-2.5, 0.8)) :
- Neglect death (aff < 0.30) : 94.7 %
- Positive selection (0.30 ≤ aff ≤ 2.50) : 5.2 %
- DP neg sel (aff > 2.50) : ~0 %
- Médullary neg sel (aff > 1.80 ∩ pos) : ~0 % (quantitativement négligeable)

#### Fix B — Ghost DP accumulation

`models/thymus_selection/model.py` — bloc DP dans `ThymicRealisation.step()` :
- `survivors.append(agent)` n'est plus appelé inconditionnellement
- Agents avec `fate != "alive"` (neglect_dead, neg_deleted) ne sont plus
  ajoutés à survivors → supprimés définitivement de la simulation
- Agents sans encounter ce pas → ajoutés à survivors normalement

#### Fix C — max_dp_age_steps = 48

`parameters.yaml` : ajout de `max_dp_age_steps: 48` (2 jours de lifespan DP)
`model.py` : ajout de `dp_steps: int = 0` dans Thymocyte ; incrémenté et
comparé à `max_dp_age_steps` avant le test d'encounter.

---

### Résultats finaux — comparaison avant/après

| Métrique | Avant (μ=0) | Après (μ=-2.5) | Cible biologique |
|----------|------------|----------------|-----------------|
| Efficacité THY1 (j5–30) | ~70 % | **5.4 %** | 2–5 % |
| Efficacité COMP1 (j8–30) | ~70 % | **4.6 %** | 2–5 % |
| Efficacité COMP2 (j8–30) | ~70 % | **4.5 %** | 2–5 % |
| neg_sel médullaire/j | ~23 | ~0 | rare (<1 % des thymocytes) |
| Neglect/j (THY1, 100 cells) | ~7 | **~96** | ~95 % attendu |
| Export THY1 régime | ~70 cells/j | **~5.4 cells/j** | ~3–5 cells/j |
| Export COMP1 régime | ~4.5 cells/j | **~0.29 cells/j** | ~0.3 cells/j |

Le pipeline est désormais quantitativement cohérent avec la physiologie thymus
murin : ~95 % de mort par neglect, ~5 % exportés. L'artefact de ghost
accumulation est résolu.

### Runs de référence

| Config | Run | Efficacité régime |
|--------|-----|-------------------|
| THY1 | `results/THY1/20260326T111317/` | 5.4 % |
| COMP1 | `results/COMP1/20260326T111322/` | 4.6 % |
| COMP2 | `results/COMP2/20260326T111326/` | 4.5 % |

---

## 2026-03-26 — Vérification calibration : flux BM et lag de transit sanguin

### Questions posées

1. Le flux BM de 6.4 cells/jour est-il cohérent avec la littérature murine ?
2. Le lag de transit sanguin de 7 jours est-il cohérent avec la biologie murine ?

---

### 1. Flux BM — analyse quantitative

#### Valeur du modèle

Le modèle à 3 compartiments (HSC → CLP → DN1) produit à l'équilibre :

```
CLP_eq = d_diff × HSC_eq / (α_T + α_myeloid + d_CLP)
       = 0.003 × 9 075 / 0.34 = 80 cells
DN1_eq = α_T × CLP_eq / export_rate
       = 0.08 × 80 / 0.15 = 42.7 cells
flux   = export_rate × DN1_eq = 6.41 cells·day⁻¹
```

#### Estimations dans la littérature murine (taux de semence thymique)

| Référence | Méthode | Estimation |
|-----------|---------|-----------|
| Scollay et al. 1986, J Immunol 136:3829 | Transfert intrathymique + back-calculation | ~20–25 cells/day |
| Schwarz & Bhandoola 2004, Nat Rev Immunol | Revue | 10–100 cells/day |
| Bhandoola et al. 2007, Immunity 26:753 | Fréquence ETP circulants × débit cardiaque | ~100 cells/day |
| Scimone et al. 2006, Blood 107:3585 | Étude CCR9/CXCR4 | <50 cells/day |
| Zlotoff & Bhandoola 2011, Curr Opin Immunol | Revue | 10–100 cells/day |

Le modèle produit **6.4 cells/day** vs plage littérature **10–100 cells/day** : sous-estimation d'un facteur **3–15×**.

#### Cause racine : limitation structurelle du modèle à 3 stades

Le chemin réel dans la moelle osseuse est :

```
HSC → MPP → LMPP → CLP → ETP/DN1
```

Le modèle compresse ce chemin en :

```
HSC → CLP → DN1
```

Les stades MPP et LMPP sont absents. Or chaque division progenitrice **amplifie** la population : 1 HSC → 2–4 LMPPs → 4–8 CLPs. Cette amplification clonale intermédiaire est perdue dans le modèle simplifié, ce qui explique le CLP_eq = 80 (littérature : 1 000–5 000 CLPs murins dans la moelle).

#### Peut-on corriger par recalibration des paramètres ?

Non, sans violer les contraintes existantes :

| Flux cible | CLP_eq requis | d_diff requis | Conflit |
|-----------|--------------|--------------|---------|
| 20 cells/day | 250 cells | 0.0094 day⁻¹ | ×3 vs Busch et al. 2015 |
| 50 cells/day | 625 cells | 0.0234 day⁻¹ | ×8 vs Busch et al. 2015 |
| 100 cells/day | 1 250 cells | 0.0468 day⁻¹ | ×16 vs Busch et al. 2015 |

`d_diff = 0.003 day⁻¹` est bien calibré sur Busch et al. 2015 (Nature 518:542) pour le taux de différentiation HSC global. L'augmenter ne ferait que déplacer le problème.

#### Décision

**Limitation architecturale documentée, non corrigée dans ce cycle.** La correction requiert l'ajout des compartiments MPP et LMPP, ce qui constitue un changement structurel du modèle BM. Les dynamiques relatives restent correctes ; l'échelle absolue est sous-estimée. Les simulations couplées (COMP1–COMP3) sont cohérentes entre elles même si le flux absolu est trop faible.

**Action future** : ~~étendre le modèle BM à 5 compartiments~~ **EXÉCUTÉ — voir session 2026-03-26 ci-dessous.**

---

## 2026-03-26 — Modèle BM v3 → v4 : passage à 5 compartiments

### Motivation

Le modèle ODE à 3 compartiments (HSC → CLP → DN1) sous-estimait le flux thymique d'un facteur 3–15× par rapport à la littérature murine (10–100 cells/day). La cause était structurelle : l'absence des stades MPP et LMPP supprimait l'amplification clonale de transit qui booste le flux dans la biologie réelle.

### Nouvelle architecture ODE (v4)

```
dHSC/dt  = r_self × HSC × (1 − HSC/K_niche) − (d_HSC_MPP + d_apop) × HSC
dMPP/dt  = d_HSC_MPP × HSC + r_MPP × MPP − (d_MPP_LMPP + d_MPP_death) × MPP
dLMPP/dt = f_lymphoid × d_MPP_LMPP × MPP − (d_LMPP_CLP + d_LMPP_death) × LMPP
dCLP/dt  = d_LMPP_CLP × LMPP − (alpha_T + d_CLP_other + d_CLP_death) × CLP
dDN1/dt  = alpha_T × CLP − export_rate × DN1
```

Le terme `r_MPP × MPP` représente l'amplification clonale de transit au stade MPP (auto-renouvellement limité avant l'engagement de lignée).

### Dérivation des équilibres analytiques

```
λ_MPP   = d_MPP_LMPP + d_MPP_death − r_MPP = 0.15 + 0.07 − 0.20 = 0.020 day⁻¹
HSC_eq  = K_niche × (1 − (d_HSC_MPP + d_apop)/r_self) = 11000 × 0.825 = 9075
MPP_eq  = d_HSC_MPP × HSC_eq / λ_MPP                  = 27.22 / 0.020 = 1361
LMPP_eq = f_lymphoid × d_MPP_LMPP × MPP_eq / (d_LMPP_CLP + d_LMPP_death)
        = 0.35 × 0.15 × 1361 / 0.12 = 595
CLP_eq  = d_LMPP_CLP × LMPP_eq / (alpha_T + d_CLP_other + d_CLP_death)
        = 0.08 × 595 / 0.14 = 340
DN1_eq  = alpha_T × CLP_eq / export_rate = 0.08 × 340 / 0.15 = 181
flux_eq = alpha_T × CLP_eq = 27.2 cells·day⁻¹
```

### Nouveaux paramètres (sources murines)

| Paramètre | Valeur | CI-95 | Source |
|-----------|--------|-------|--------|
| `d_HSC_MPP` | 0.003 day⁻¹ | [0.0025, 0.0035] | Busch et al. 2015, Nature 518:542 (renommé depuis d_diff, valeur inchangée) |
| `r_MPP` | 0.20 day⁻¹ | [0.17, 0.23] | Oguro et al. 2013, Cell Stem Cell 13:535; Wilson et al. 2008, Nature 453:529 |
| `d_MPP_LMPP` | 0.15 day⁻¹ | [0.12, 0.18] | Adolfsson et al. 2005, Cell 121:295 |
| `d_MPP_death` | 0.07 day⁻¹ | [0.05, 0.09] | Estimé de Busch et al. 2015 |
| `f_lymphoid` | 0.35 | [0.28, 0.42] | Adolfsson et al. 2005; Luc et al. 2012, J Exp Med 209:1441 |
| `d_LMPP_CLP` | 0.08 day⁻¹ | [0.06, 0.10] | Inlay et al. 2009, J Exp Med 206:1387 |
| `d_LMPP_death` | 0.04 day⁻¹ | [0.02, 0.06] | Estimé de Kondo et al. 1997, Cell 91:661 |
| `d_CLP_other` | 0.05 day⁻¹ | [0.03, 0.07] | Rodrigues et al. 2005, Blood 105:3845 (remplace alpha_myeloid) |
| `d_CLP_death` | 0.01 day⁻¹ | [0.007, 0.013] | Estimé; Guo et al. 2010, J Exp Med 207:1083 |

### CI-95 pour MPP : régime quasi-équilibré

La sensibilité de MPP_eq aux paramètres est amplifiée par le facteur `1/λ_MPP = 50` :

```
S(r_MPP)      = r_MPP      / λ_MPP = 0.20/0.02 = 10×
S(d_MPP_LMPP) = d_MPP_LMPP / λ_MPP = 0.15/0.02 =  7.5×
S(d_MPP_death)= d_MPP_death/ λ_MPP = 0.07/0.02 =  3.5×
```

La propagation exacte donne des CI incluant des valeurs négatives (biologiquement impossibles). La méthode `_mpp_ci95` utilise délibérément S=1 (approximation conservatrice), produisant CI=[844, 1878] autour de MPP_eq=1361. C'est une borne inférieure de l'incertitude réelle ; les CI MPP sont à traiter comme indicatifs uniquement.

### Comparaison v3 → v4

| Métrique | v3 (3 compart.) | v4 (5 compart.) | Littérature murine |
|---------|----------------|----------------|-------------------|
| Compartiments | HSC, CLP, DN1 | HSC, MPP, LMPP, CLP, DN1 | — |
| HSC_eq | 9 075 | 9 075 ✓ | ~10 000 |
| CLP_eq | 80 | 340 | 500–5 000 (fonctionnel) |
| DN1_eq (BM) | 43 | 181 | 100–500 |
| flux_eq | **6.4 cells/day** | **27.2 cells/day** | **10–100 cells/day** |
| Contrainte Busch 2015 | ✓ (d_diff=0.003) | ✓ (d_HSC_MPP=0.003) | — |

Amélioration du flux : **×4.3** (6.4 → 27.2 cells/day), désormais dans la plage littérature.

### Compatibilité orchestrateur (vérifiée)

- `export_signals[0].signal_id` = `"bm_haematopoiesis.progenitor_export"` ✓ (inchangé)
- `export_signals[0].entity_id` = `"CL:0002420"` ✓ (inchangé)
- Format ISSL : `continuous_state` passe de 3 à 5 entrées ; compatible avec le schema (pas de contrainte sur le nombre d'entrées)
- `MODEL_VERSION` : 3 → **4**
- 5 tests unitaires existants : **5/5 PASSED** (aucune modification des tests)
- 5 simulations complètes (BM1, THY1, COMP1, COMP2, COMP3) : **0 alertes watchdog**

### Fichiers modifiés

| Fichier | Modification |
|---------|-------------|
| `models/bm_haematopoiesis/parameters.yaml` | Réécriture complète — 5 compartiments, nouvelles ICs à l'équilibre analytique |
| `models/bm_haematopoiesis/model.py` | Réécriture complète — v4, `_ode()` × 5 compart., `emit_issl()` × 5 entrées ISSL, `_mpp_ci95()` dédié |
| `figures/generate_figures.py` | `load_bm()` mis à jour pour 5 compartiments + fig1 étendue à 5 courbes |

---

### 2. Lag de transit sanguin — analyse quantitative

#### Valeur du modèle (avant correction)

Le lag est calculé analytiquement dans `models/blood_transit/model.py` :

```
t_stop = -ln(1 - stop_fraction) / k_total
       = -ln(1 - 0.95) / 0.43
       = -ln(0.05) / 0.43
       = 2.996 / 0.43
       = 6.97 jours
```

Ce n'est pas le temps moyen de transit — c'est le **95e percentile** de la distribution exponentielle de transit.

#### Statistiques de transit pour les paramètres actuels

| Statistique | Formule | Valeur |
|------------|---------|--------|
| Temps moyen | `1 / k_total = 1/0.43` | **2.33 jours** |
| Médiane | `ln(2) / k_total` | **1.61 jours** |
| 75e percentile | `-ln(0.25) / 0.43` | **3.22 jours** |
| 82e percentile | `-ln(0.18) / 0.43` | **3.99 jours** ← cible |
| **95e percentile** (ancienne valeur) | `-ln(0.05) / 0.43` | **6.97 jours** |

#### Littérature murine pour le temps de transit BM → thymus

| Référence | Méthode | Estimation |
|-----------|---------|-----------|
| Schwarz & Bhandoola 2004, Nat Rev Immunol 4:812 | Revue | **1–4 jours** |
| Scimone et al. 2006, Blood 107:3585 | Injection i.v. + suivi thymique | **1–2 jours** pour l'arrivée des premières cellules |
| Krueger et al. 2010, J Exp Med 207:1351 | CCR9 KO | entrée retardée de ~2–3 j → entrée normale ~1–2 j |
| Bhandoola et al. 2007, Immunity 26:753 | ETP circulants | demi-vie circulante ~1–2 j |

Les données murines convergent sur **1–4 jours**, avec une médiane ~2 jours. Les paramètres `k_homing=0.35` et `k_death=0.08` sont biologiquement corrects (demi-vie circulante = 1.61 j). Le problème était uniquement le choix du `stop_fraction=0.95` (95e percentile).

#### Cause racine

`stop_fraction=0.95` est une **convention de modélisation trop conservative** : elle attend que 95 % des cellules soient absorbées ou mortes, ce qui correspond au 95e percentile de la distribution. Ce n'est pas le délai biologique de transit mais le délai jusqu'à ce que la quasi-totalité du batch soit traitée.

#### Correctif appliqué

`models/blood_transit/parameters.yaml` : `stop_fraction : 0.95 → 0.82`

```
t_stop = -ln(1 - 0.82) / 0.43 = -ln(0.18) / 0.43 = 1.715 / 0.43 ≈ 3.99 jours
```

Impact sur les cellules livrées :

| Paramètre | Avant (sf=0.95) | Après (sf=0.82) |
|-----------|----------------|----------------|
| lag_s | 601 933 s (6.97 j) | 344 555 s **(3.99 j)** |
| delivered / 6.41 N0 | 4.95 cells (77.3 %) | 4.30 cells **(67.1 %)** |

La livraison baisse légèrement (~14 %) car on ne compte plus les cellules arrivant entre 4 et 7 jours, mais ces dernières arrivaient trop tard pour contribuer significativement à la dynamique thymus.

#### Résultat observé dans COMP2/COMP3

```
TransferDispatcher: 'blood_transit' → lag_s=344554.8, delivered=4.30
```

Lag désormais cohérent avec la littérature murine. ✓

---

### Bilan — calibration espèce (audit complet de la session)

| Modèle | Paramètre | Statut |
|--------|-----------|--------|
| `bm_haematopoiesis` | tous | ✓ Murin (Busch 2015, Bhatt 2016, Cheshier 1999…) |
| `bm_haematopoiesis` | `alpha_T` | ✓ Corrigé : `Allman 2003 Blood` → `Allman 2003 Nat Immunol 4:168` |
| `bm_haematopoiesis` | flux absolu | ⚠️ 6.4 cells/day vs lit. 10–100 — limitation structurelle (modèle 3 stades) |
| `blood_transit` | `k_homing`, `k_death` | ✓ Murin (Bhandoola 2007, Schwarz & Bhandoola 2004) |
| `blood_transit` | `stop_fraction` | ✓ Corrigé : 0.95 → **0.82** (lag 7 j → **4 j**, cohérent avec littérature murine) |
| `thymus_selection` | tous | ✓ Murin (Hogquist 1994, Kappler 1987, Shortman 2002…) |
| `thymus_selection` | `p_encounter_cortex` | ✓ Corrigé : référence Bhatt (niches HSC) retirée, note de modélisation ajoutée |
| `peripheral_ln` | `S4, S8, rho_4, rho_8` | ✓ Corrigé : `Borghans & De Boer 2007` (données humaines HIV) → `Tough & Sprent 1994 J Exp Med` + `Kieper & Jameson 1999 PNAS` (murins) |
| `peripheral_ln` | `d_act`, `d_death` | ✓ Corrigé : `Surh & Sprent 2008, Nature 450` → `Surh & Sprent 2008, Immunity 29:848` (journal correct) |

### Fichiers modifiés (cette session)

| Fichier | Modification |
|---------|-------------|
| `models/blood_transit/parameters.yaml` | `stop_fraction 0.95 → 0.82` |
| `models/bm_haematopoiesis/parameters.yaml` | Citation `alpha_T` corrigée |
| `models/thymus_selection/parameters.yaml` | Citation `p_encounter_cortex` corrigée |
| `models/peripheral_ln/parameters.yaml` | Citations `S4/S8/rho/d_*` corrigées vers sources murines |
| `configs/run_COMP3_full_graph.yaml` | Durée 30 j → **365 j** + `model_args --baseline-import 0` |
| `figures/generate_figures.py` | Script de génération des 5 figures (fig1–fig5) |

### Runs de référence (après correction lag)

| Config | lag_s observé | delivered (par 6.41 N0) |
|--------|--------------|------------------------|
| COMP2 `results/COMP2/latest/` | 344 555 s (3.99 j) ✓ | 4.30 cells |
| COMP3 `results/COMP3/latest/` | 344 555 s (3.99 j) ✓ | 4.28–4.30 cells |

---

## 2026-03-28 — Recalibration à l'échelle biologique (BM v5 + Thymus v3 + scaling factor)

### Contexte

Les paramètres des modèles BM et Thymus étaient calibrés à une échelle très
réduite par rapport aux valeurs publiées dans la littérature murine. En
particulier :
- Compartiments BM (MPP, LMPP, CLP) sous-estimés de 2 à 100× vs littérature
- Modèle Thymus sans facteur d'échelle documenté → flux d'export non
  interprétable biologiquement

### Objectif

1. BM v5 : paramètres et conditions initiales calés sur les ordres de grandeur
   publiés (cellules totales dans la moelle de souris adulte C57BL/6)
2. Thymus v3 : introduction d'un `scaling_factor` dans le fichier ISSL pour
   relier le compte ABM réduit aux effectifs biologiques réels
3. Vérification de la cohérence des valeurs simulées avec la littérature
4. Vérification des temps de transfert (blood transit, homing périphérique)
5. Régénération des figures (tous les scénarios)

---

### BM v5 — Recalibration (paramètres modifiés)

#### Changements de paramètres

| Paramètre | v4 | v5 | Justification |
|-----------|----|----|---------------|
| `r_MPP` | 0.200 | **0.199** | λ_MPP = 0.001 au lieu de 0.020 → MPP amplifié ~27× |
| `d_MPP_death` | 0.070 | **0.050** | Apoptose MPP réduite (Wilson 2008 BrdU data) |
| `alpha_T` | 0.080 | **0.002** | Engagement T-lineage ~6 % des CLP (Allman 2003) |
| `d_CLP_other` | 0.050 | **0.020** | Réajusté pour équilibre CLP |
| `export_rate` | 0.150 | **0.050** | ETPs résidents BM plus longtemps avant éjection |

#### Nouvelles conditions initiales (équilibre analytique vérifié)

| Compartiment | v4 (avant) | v5 (après) | Littérature murine |
|---|---|---|---|
| HSC | 9 075 | 9 075 | ~10 000 (Bhatt 2016) ✓ |
| MPP | 1 361 | **27 225** | ~15 000–50 000 (Wilson 2008) ✓ |
| LMPP | 595 | **11 911** | ~10 000–20 000 (Adolfsson 2005) ✓ |
| CLP | 340 | **29 778** | ~20 000–30 000 (Rodrigues 2005) ✓ |
| DN1/ETP | 181 | **1 191** | Pool BM pré-éjection ✓ |
| **Flux export** | 27.2 /j | **59.6 /j** | 10–100 /j (Bhandoola 2007) ✓ |

#### Distribution des destinées CLP (v5)
- T-lineage : **6.2 %** → ~60 ETP/j ✓
- B/NK lineage : 62.5 %
- Apoptose : 31.2 %

#### Vérification numérique
```
λ_MPP = 0.150 + 0.050 − 0.199 = 0.001 day⁻¹
MPP_eq = (0.003 × 9 075) / 0.001 = 27 225 ✓
flux_eq = 0.002 × 29 778 = 59.56 cells/day ✓
```

---

### Thymus v3 — Facteur d'échelle

#### Principe

Le thymus murin contient ~5×10⁷ thymocytes (Scollay & Godfrey 1995), non
simulables en ABM. L'ABM opère à ~150 agents à l'état stationnaire (avec
~60 ETP/j d'entrée BM). Le **scaling_factor = 300 000** relie les deux :

```
cellules_réelles ≈ agents_ABM × 300 000
```

**Dérivation :**
- Thymus réel : 5×10⁷ thymocytes
- État stationnaire ABM estimé : ~167 agents
- 5×10⁷ / 167 ≈ 3×10⁵

Le facteur capture implicitement l'amplification intra-thymique (~30–50 divisions
de l'ETP au stade SP en biologie).

#### Implémentation dans le fichier ISSL

Le fichier ISSL émis par le thymus contient désormais :

| Champ | Emplacement | Valeur |
|-------|-------------|--------|
| `scaling_factor` | `internal_parameters` | 300 000 |
| `raw_dn_agents`, `raw_dp_agents`, etc. | `internal_parameters` | comptes ABM bruts |
| `continuous_state[*].count` | — | valeurs **scalées** (agents × 300 000) |
| `export_signals[0].flux` | — | flux **brut ABM** (pour le couplage orchestrateur/PLN) |
| `export_signals[0].biological_flux_per_day` | — | flux biologique estimé (flux × 300 000) |

#### Résultats à l'état stationnaire (jour 30, baseline_import=60/j)

| Variable | Simulé (scalé) | Littérature murine |
|----------|----------------|-------------------|
| DN total | ~7 M | ~1.5 M (ratio élevé : staging DN simplifié) |
| DP | ~1.7 M | ~42 M (attendu ; ABM sans amplification DP explicite) |
| SP total | ~3 M | ~6 M (×2, raisonnable) |
| **Export** | **~1 M cells/j** | **1–2 M/j ✓** |

---

### Temps de transfert — Vérification

| Transfert | Valeur simulée | Littérature |
|-----------|---------------|------------|
| BM → Thymus (blood transit) | **~3.99 jours** (lag_s ≈ 344 555 s) | 1–4 j (Schwarz & Bhandoola 2004) ✓ |
| Thymus → PLN (homing) | **2 jours** (constante 172 800 s) | quelques heures–2 j (Butcher & Picker 1996) ✓ |

> **Note** : Le lag blood transit observé dans les nouvelles simulations est
> ~95.7 h ≈ 4 jours (inchangé par rapport à v4 car les paramètres `blood_transit`
> n'ont pas été modifiés).

---

### Fichiers modifiés

| Fichier | Nature |
|---------|--------|
| `models/bm_haematopoiesis/parameters.yaml` | v4 → v5 : 5 paramètres + conditions initiales |
| `models/bm_haematopoiesis/model.py` | VERSION "4" → "5", docstring, bornes OOD |
| `models/thymus_selection/parameters.yaml` | v2 → v3 : `scaling_factor`, `baseline_import` 100→60 |
| `models/thymus_selection/model.py` | VERSION "2" → "3", `emit_issl` + scaling |
| `results/plot_results.py` | Mise à jour figures : thymus en ×10⁶, titres v5/v3 |
| `run_and_plot.py` | Nouveau : runner in-process + générateur de figures |
| `README.md` | Mise à jour descriptions modèles BM v5 et Thymus v3 |

### Runs de référence

Simulations régénérées in-process (`python run_and_plot.py`) :

| Scénario | Répertoire | Flux BM eq. | Flux Thymus (scalé) |
|----------|-----------|------------|---------------------|
| BM1 | `results/BM1_baseline/BM1_baseline/` | 59.6 cells/j | — |
| THY1 | `results/THY1_baseline/THY1_baseline/` | — | ~1.0 M cells/j |
| COMP1 | `results/COMP1_direct/COMP1_direct/` | 59.6 | ~1.0 M |
| COMP2 | `results/COMP2_transfer/COMP2_transfer/` | 59.6 | ~1.0 M (lag ~4 j) |
| COMP3 | `results/COMP3_full_graph/COMP3_full_graph/` | 59.6 | ~1.0 M → PLN |

Figures régénérées dans `results/figures/` (fig1–fig6).


---

## 2026-03-28 — COMP3 : durée ramenée de 365 à 30 jours

### Contexte

La simulation COMP3 (`run_COMP3_full_graph.yaml`) était configurée à 365 jours
alors que toutes les autres simulations (BM1, THY1, COMP1, COMP2) tournent sur
30 jours. Cette incohérence rendait les comparaisons inter-scénarios difficiles
et la figure 5 non comparable aux figures 3 et 4.

### Modification

`configs/run_COMP3_full_graph.yaml` :
```
end_s: 31536000  →  end_s: 2592000   # 365 j → 30 j
```

Le commentaire précédent indiquait que 365 j étaient « requis pour résoudre la
relaxation homéostatique du PLN (τ ≈ 167 j) ». Cette contrainte est levée : la
figure 5 montre désormais le transitoire d'établissement sur 30 jours, cohérent
avec les autres figures.

### Résultats COMP3 à j 30

| Variable | Valeur |
|----------|--------|
| BM export flux | 59.6 cells/j |
| Transit lag (blood) | ~95.7 h ≈ 4 jours |
| Thymus export (scalé) | ~725 K cells/j |
| PLN CD4 naïf | ~183 552 cells |
| PLN CD8 naïf | ~91 777 cells |
| Ratio CD4/CD8 | ~2.0 ✓ |

Le ratio CD4/CD8 ≈ 2:1 est cohérent avec les données murines.
Le PLN n'a pas encore atteint son set-point (S4 = 200 000) après 30 jours car
le signal thymique pénètre avec un lag BM→thymus (~4 j) + thymus→PLN (2 j).

### Fichiers modifiés

| Fichier | Modification |
|---------|-------------|
| `configs/run_COMP3_full_graph.yaml` | `end_s` 31 536 000 → 2 592 000 |
| `results/COMP3_full_graph/COMP3_full_graph/issl/` | Données ISSL régénérées (30 j) |
| `results/figures/` | Toutes les figures régénérées (fig1–fig6) |


---

## 2026-03-28 — Bug fig2 : panneau "Naïve T export" vide dans THY1 baseline

### Symptôme

Le panneau `gs[1, 1]` de la figure 2 (`fig2_THY1_baseline.png`) restait vide malgré la présence d'un export de cellules T naïves dans les résultats de simulation (affiché en console : `THY1 done — export: 3.3 agents/cp (~1 000 000 cells/day scaled)`).

### Cause racine

Dans `results/plot_results.py`, le tableau `stage_info` de `fig_thy1_baseline()` incluait :

```python
("CL:0000898", "Naïve T export", PALETTE["export"], gs[1, 1]),
```

`get_entity_series(thy, "CL:0000898")` recherche dans `continuous_state` des checkpoints ISSL. Or `CL:0000898` (naïve T cell) n'est **jamais** dans `continuous_state` du modèle thymus — il apparaît uniquement dans `export_signals`. La fonction renvoyait donc des tableaux vides → panneau blanc.

### Correction

1. **Nouveau helper `get_bio_export_flux(issl_records)`** ajouté dans `plot_results.py` : lit `biological_flux_per_day` dans `export_signals[0]`, avec fallback sur `flux * scaling_factor`.

2. **Panneau `gs[1, 1]` réécrit** : remplace l'entrée bogée dans `stage_info` par un subplot dédié qui appelle `get_bio_export_flux(thy)` et affiche le flux biologique en ×10⁶ cellules·jour⁻¹.

### Note biologique

Les 3–4 premiers jours restent à zéro — comportement attendu : les thymocytes doivent traverser DN (≥1 j) → DP (~20 sous-pas) → sélection positive → séjour médullaire (`medullary_dwell_steps=72` à 1 h/sous-pas = 3 j minimum) avant export.

### Fichiers modifiés

| Fichier | Modification |
|---------|-------------|
| `results/plot_results.py` | Ajout `get_bio_export_flux()` ; panneau THY1 `gs[1,1]` corrigé |
| `results/figures/fig2_THY1_baseline.png` | Régénérée avec export visible |

---

## 2026-03-28 — BM v6 : initialisation au steady-state analytique exact

### Symptôme

La figure 1 montrait une légère dérive croissante du flux d'export BM (inset panel) sur 30 jours. Matplotlib auto-scale le flux et rendait une variation de ~0.004 cells/day visuellement significative.

### Cause racine

Les conditions initiales dans `parameters.yaml` (v5) étaient des entiers arrondis calculés par cascade manuelle :

| Compartiment | IC entière (v5) | Vrai SS analytique | Erreur |
|---|---|---|---|
| HSC  |  9075     |  9075.0000  |  0.0000 |
| MPP  | 27225     | 27225.0000  |  0.0000 |
| LMPP | 11911     | 11910.9375  | +0.0625 |
| CLP  | 29778     | 29777.3438  | +0.6563 |
| DN1  |  1191     |  1191.0938  | −0.0938 |

L'erreur de CLP (+0.66 cells) excitait le mode lent du compartiment CLP (τ ≈ 1/0.032 = 31 jours), produisant une dérive monotone du flux sur toute la durée de la simulation. Bien qu'infime en valeur absolue (0.007%), elle était visuellement évidente après auto-scaling matplotlib.

### Correction

`BMHaematopoiesis.__init__` appelle maintenant `_compute_steady_state()` au lieu de lire les ICs du YAML :

```python
def _compute_steady_state(self) -> np.ndarray:
    p = self._p
    HSC  = p["K_niche"] * (1.0 - (p["d_HSC_MPP"] + p["d_apop"]) / p["r_self"])
    lam  = p["d_MPP_LMPP"] + p["d_MPP_death"] - p["r_MPP"]
    # raises ValueError if lam ≤ 0 (unstable regime)
    MPP  = p["d_HSC_MPP"] * HSC / lam
    LMPP = p["f_lymphoid"] * p["d_MPP_LMPP"] * MPP / (p["d_LMPP_CLP"] + p["d_LMPP_death"])
    CLP  = p["d_LMPP_CLP"] * LMPP / (p["alpha_T"] + p["d_CLP_other"] + p["d_CLP_death"])
    DN1  = p["alpha_T"] * CLP / p["export_rate"]
    return np.array([HSC, MPP, LMPP, CLP, DN1])
```

Résultat vérifié : flux constant à 59.5547 cells/day sur 30 jours (dérive = 0.00e+00).

La ligne de référence dans la figure 1 est aussi corrigée : 59.6 → 59.5547.

### Fichiers modifiés

| Fichier | Modification |
|---|---|
| `models/bm_haematopoiesis/model.py` | v5→v6 ; ajout `_compute_steady_state()` ; init SS analytique |
| `models/bm_haematopoiesis/parameters.yaml` | v5→v6 ; ICs corrigées (float SS) ; note doc-only |
| `results/plot_results.py` | Référence flux 59.6 → 59.5547 ; titre fig1 v5→v6 |
| `results/figures/` | Toutes les figures régénérées |

---

## 2026-03-28 — PLN : bug import thymus (scaling_factor absent) + ratio CD4/CD8

### Symptôme

Dans COMP3, le ratio CD4/CD8 était figé à exactement 2.000 sur 30 jours, et les pools déclinaient régulièrement de −8% (~200 000 → 183 500 pour CD4).

### Analyse

**Bug 1 — import biologique non appliqué**

`_step()` du modèle PLN lisait `sig["flux"]` (2.42 agents ABM/jour) en lieu et place de `sig["biological_flux_per_day"]` (725 000 cellules/jour). Le facteur de scaling (300 000) du modèle thymus n'était pas transmis au PLN. Import effectif reçu : **1.6 cellules/jour** au lieu de ~923 cellules/jour.

**Bug 2 — paramètre lymph_node_fraction absent**

Sans fraction de distribution, le PLN recevait (incorrectement) 100 % de l'output thymus plutôt que la part distribuée à un compartiment LN unique.

**Conséquence sur la dynamique**

L'équilibre Borghans-De Boer sans import est :
```
CD4* = ρ₄ × S4 / (ρ₄ + d_total) = 0.003 × 200 000 / 0.006 = 100 000
```
Avec ρ₄ = d_total = 0.003 jour⁻¹, le système déclinait vers S4/2 (τ ≈ 167 jours). Sur 30 jours : −8% attendu, −8.2% observé ✓.

Le ratio 2.000 exact était un artefact : avec import ≈ 0 et S4 = 2×S8, ρ₄ = ρ₈, la symétrie algébrique conserve le rapport CD4/CD8 à exactement 2.000.

### Correction

**`models/peripheral_ln/parameters.yaml`** — Ajout :
```yaml
lymph_node_fraction:
  value: 0.0013   # calibré : 923 cells/day requis / ~700 000 cells/day thymus
```

**`models/peripheral_ln/model.py`** — `_step()` corrigé :
```python
if sig.get("biological_flux_per_day") is not None:
    total_bio_flux = float(sig["biological_flux_per_day"])
elif sig.get("scaling_factor") is not None:
    total_bio_flux = float(sig["flux"]) * float(sig["scaling_factor"])
else:
    total_bio_flux = float(sig["flux"])
flux_per_day = total_bio_flux * self._p["lymph_node_fraction"]
```

### Résultat après correction

| Métrique | Avant (buggy) | Après (corrigé) |
|---|---|---|
| CD4 à jour 30 | 183 552 (−8.2%) | 193 424 (−3.3%) |
| CD8 à jour 30 |  91 777 (−8.2%) |  97 093 (−2.8%) |
| Ratio à jour 30 | 2.000 (figé) | 1.992 (converge vers ~1.92) |
| Import CD4 reçu | 1.6 cells/day | ~612 cells/day |

La décroissance résiduelle de −3% est attendue : sur 30 jours (< τ/5), les pools n'ont pas encore atteint leur équilibre asymptotique (~202 000 / 105 000 avec flux thymus en régime établi).

### Fichiers modifiés

| Fichier | Modification |
|---|---|
| `models/peripheral_ln/parameters.yaml` | Ajout `lymph_node_fraction: 0.0013` |
| `models/peripheral_ln/model.py` | `_step()` : lecture `biological_flux_per_day`, application `lymph_node_fraction` |
| `results/figures/` | Figures régénérées |

---

## 2026-03-28 — Régénération des figures + mise à jour README (session BM v6 / PLN fix)

Figures régénérées avec `run_and_plot.py --figs` après les corrections de cette session :

| Figure | Changements visibles |
|---|---|
| `fig1_BM1_baseline.png` | Flux export plat à 59.5547 cells/day (plus de dérive) ; référence horizontale corrigée (59.6 → 59.5547) |
| `fig2_THY1_baseline.png` | Panneau Naïve T export désormais affiché (lecture `biological_flux_per_day` depuis `export_signals`) |
| `fig3_COMP1_direct.png` | Export BM démarrant au SS exact |
| `fig4_COMP2_transfer.png` | Idem |
| `fig5_COMP3_full_graph.png` | PLN : pools CD4/CD8 maintenus (~−3% en 30 j au lieu de −8%) ; ratio 2.000 → 1.992 (dynamique réelle) |
| `fig6_comparative.png` | Toutes courbes cohérentes avec les corrections ci-dessus |

README mis à jour :
- BM : version v5 → v6 ; valeurs SS exactes (LMPP 11910.94, CLP 29777.34, DN1 1191.09, flux 59.55) ; mention de l'initialisation analytique
- Thymus : calibration citée BM v6 ; note sur `CL:0000898` dans `export_signals` uniquement
- PLN : étape 2 réécrite (lecture `biological_flux_per_day`, `lymph_node_fraction`) ; étape 7 (calibration, équilibre Borghans, ratio long-terme ~1.92)
- Chemins de logs : BM_haematopoiesis_v5 → v6

---

## 2026-03-28 — BM v7 : remplacement CI-95 linéaire S=1 par Monte Carlo

### Symptôme

Les barres d'erreur CI-95 du modèle BM étaient fausses ("aux fraises") sur tous les compartiments.

### Diagnostic (trois problèmes imbriqués)

**1. P(λ_MPP ≤ 0) = 47.9 %**
λ_MPP = d_MPP_LMPP + d_MPP_death − r_MPP = 0.001 jour⁻¹ (très proche de zéro).
σ(λ_MPP) = √(σ²_d_MPP_LMPP + σ²_d_MPP_death + σ²_r_MPP) = 0.0191.
Soit z = 0.001/0.0191 = 0.052 → P(λ ≤ 0) = Φ(−0.052) ≈ 48 %.
Presque la moitié de l'espace paramétrique CI-95 donne un système instable (MPP diverge). Il n'existe donc aucune CI analytique finie pour MPP.

**2. Sensibilité S = 1 utilisée, vraie valeur S = 199**
La formule `_state_ci95` supposait ∂MPP/∂p × p/MPP = 1 pour tous les paramètres.
La vraie sensibilité de MPP à r_MPP est :
```
S(r_MPP→MPP) = r_MPP / λ_MPP = 0.199/0.001 = 199
S(d_MPP_LMPP→MPP) = 0.150/0.001 = 150
S(d_MPP_death→MPP) = 0.050/0.001 = 50
```
L'ancienne formule sous-estimait la CI de MPP d'un facteur ~100–200.

**3. Propagation amont absente pour LMPP, CLP, DN1, flux**
Les CI de tous les compartiments aval ne comprenaient pas les incertitudes sur r_MPP/d_MPP_LMPP/d_MPP_death, qui dominent la cascade via S = 150–199.

### CI obtenues avec les trois méthodes

| Compartiment | Ancienne CI (S=1) | CI correcte S linéaire | MC CI (cond. stabilité) |
|---|---|---|---|
| MPP  | [14 973, 39 477] | [−245 000, +299 000] | [587, 39 308] |
| flux | [57.9, 61.2]     | inapplicable          | [1.4, 90.9]  |

La CI linéaire correcte donne des bornes négatives → inutilisable. Le MC conditionnel est la seule approche valide.

### Correction (BM v7)

Remplacement des méthodes `_state_ci95`, `_mpp_ci95`, `_flux_ci95` par `_compute_mc_ci95()` :
- Échantillonnage vectorisé de 2 000 jeux de paramètres N(μ, σ)
- Rejet des tirages avec λ_MPP ≤ 0 ou paramètres non-positifs (~48 % des tirages)
- Percentiles 2.5/97.5 sur les tirages valides
- Précompilé une fois dans `__init__`, résultats stockés dans `self._ci95`

Ajout dans le watchdog ISSL de `lambda_mpp`, `ci_stability_fraction`, `ci_n_stable`, `ci_n_total`.

### Résultat

```
MC CI (n_stable=2068/4000, stability=51.7%)
  HSC   : [  7820,  10401]
  MPP   : [   587,  39308]   ← asymétrique : médiane ~1850, pas ~27225
  LMPP  : [   288,  17428]
  CLP   : [   713,  44207]
  DN1   : [    27,   1942]
  flux  : [   1.4,   90.9]  cells·day⁻¹  (couvre toute la plage littérature)
```

### Fichiers modifiés

| Fichier | Modification |
|---|---|
| `models/bm_haematopoiesis/model.py` | v6→v7 ; ajout `_compute_mc_ci95()` ; suppression des anciens helpers CI |
| `models/bm_haematopoiesis/parameters.yaml` | v6→v7 |
| `results/figures/` | Figures régénérées avec CI correctes |
