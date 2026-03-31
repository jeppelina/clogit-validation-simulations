# Choice Model Validation: Monte Carlo Simulations

Replication code for simulation studies validating conditional logit (clogit) models for assortative mating research. These simulations accompany the doctoral dissertation by Jesper Lindmarker (Linköping University, Institute for Analytical Sociology).

## Overview

The conditional logit model, rooted in McFadden's (1973) random utility theory, estimates preference parameters from discrete choice data. It is used across fields: residential choice and neighbourhood sorting (Bruch & Mare, 2012), transport mode selection, and demographic studies of partner choice and assortative mating (Jepsen & Jepsen, 2002; Haandrikman & van Wissen, 2012; Gullickson, 2021).

In partner choice applications, each individual is modelled as selecting a partner from a constructed choice set drawn from the locally eligible population. The model estimates which characteristics of potential partners and of the dyad (educational similarity, shared ancestry, residential proximity) predict union formation. Identification comes from within-individual variation: comparing the chosen partner to unchosen alternatives in the same choice set. This means the model controls for all characteristics of the focal individual that do not vary across alternatives. Conditional logit can incorporate continuous covariates like residential distance, handle multiple sorting dimensions simultaneously, and yield coefficients interpretable as revealed preference parameters conditional on the opportunity structure defined by the choice set.

These properties rest on assumptions. The model treats partner choice as one-sided: each person picks independently, with no requirement that the chosen partner reciprocates. It assumes that the choice set accurately represents the relevant partner market. It assumes independence of irrelevant alternatives (IIA). When used for mediation decomposition, adding variables sequentially to assess how much of educational homogamy operates through workplace co-membership, it further assumes that the sequential structure is not confounded by unobserved factors. These assumptions are discussed in the discrete choice literature but rarely tested quantitatively in the partner choice context.

This repository tests them. We generate synthetic partner markets under known data-generating processes (DGPs), apply conditional logit as in the empirical papers, and check whether the true preference parameters are recovered. Five modules target specific threats to identification:

1. **Group size bias** (Module 1): Does clogit confound compositional effects with preferences when ancestry groups differ in size?
2. **Identification** (Module 1b): Can clogit distinguish homophily (preferring similar education) from status maximisation (preferring high education)?
3. **Spatial confounding** (Module 2): Does residential segregation bias homophily estimates when propinquity also drives partner choice?
4. **Mediation validity** (Module 3): Is the sequential decomposition of educational homogamy into preference and opportunity components reliable under unobserved confounding?
5. **Bilateral matching bias** (Module 4): How does the mismatch between two-sided matching (reality) and one-sided choice models (clogit) affect parameter estimates?

An experimental sixth module (Module 5) compares clogit to the simulation-based boundary estimation used in Paper 2 of the dissertation. See `module5/README.md` for details and caveats.

## Key Findings

**Modules 1 to 3** test specific threats under both oracle (one-sided) and bilateral (two-sided) matching:

- **Module 1**: clogit handles group imbalance well under oracle (bias < 1%). Bilateral matching adds 33 to 47% overestimation of the homophily coefficient, roughly constant across group compositions.
- **Module 1b**: Homophily and status maximisation are distinguishable under oracle matching. The absolute value in |edu_diff| breaks collinearity with alt_edu within clogit strata.
- **Module 2**: Under oracle, spatial confounding is negligible with global choice sets and a distance control (bias < 1.5%). Under bilateral matching, the homophily coefficient is overestimated by about 55%, because bilateral filtering compounds across preference dimensions.
- **Module 3**: Under oracle, the education preference coefficient is robust to unobserved confounding (bias about 4%); the workplace coefficient absorbs the confounder. Under bilateral matching, education preference bias is about 45%, but does not interact with confounding. The two biases are approximately additive. Mediation shares should be interpreted as upper bounds.

**Module 4** provides the systematic bilateral analysis:

- A calibrated threshold sweep compares homophily and status maximisation DGPs at equivalent acceptance rates (10 to 80%).
- Homophily overestimation varies from 16% (40% acceptance) to 69% (10% acceptance).
- Status maximisation is unrecoverable under bilateral matching (66 to 96% attenuation).
- The bilateral bias converges after about 20 encounter rounds and is stable across group compositions.

**Cross-module finding**: the bilateral bias is approximately additive with other threats. Confounding, spatial structure, and group composition do not interact strongly with the bilateral mechanism. The qualitative conclusions from oracle-based analyses survive; only the quantitative magnitudes need the bilateral correction.

## Simulation Design

### Shared Settings

| Parameter | Value |
|-----------|-------|
| Population size (N) | 1,000 |
| Choice set size (J) | 30 |
| Monte Carlo repetitions (R) | 200 |
| Estimator | `survival::clogit()` with `strata(id)` |

### Data-Generating Processes

The primary DGP across all modules is **homophily**:

```
U(i,j) = -1.5 × |edu_i - edu_j| + 0.8 × same_group(i,j) + ε
```

where ε ~ Gumbel(0, 1). Education is drawn from N(12, 2). This matches the empirical specification in the dissertation papers.

Two alternative DGPs appear in Modules 1b and 4:

- **Status maximisation**: U = 0.5 × alt_edu + 0.8 × same_group + ε
- **Mixed**: U = -1.5 × |edu_diff| + 0.5 × alt_edu + 0.8 × same_group + ε

### Matching Mechanisms

Modules 1 to 3 use **oracle matching** (one-sided choice) as the baseline, plus **bilateral matching** (calibrated at 20% one-sided acceptance rate) to test whether each threat interacts with two-sided matching.

Module 4 compares three mechanisms:

| Mechanism | Description |
|-----------|-------------|
| Oracle (one-sided) | Each agent chooses their highest-utility option from the full population. This is what clogit assumes. |
| Bilateral (random encounters) | Agents meet in random encounter rounds; both must accept (U > threshold). The threshold is calibrated per DGP to achieve a target acceptance rate. |
| Deferred acceptance | Gale-Shapley algorithm: proposals, tentative acceptance, displacement. Included as a robustness check. |

The calibrated threshold approach matters because raw thresholds are meaningless across DGPs: utility distributions differ. Calibrating thresholds to achieve the same one-sided acceptance rate equalizes market selectivity for comparison.

## File Structure

```
├── README.md
├── utils.R                       Shared functions: agent generation, matching
│                                 mechanisms, choice sets, estimation, diagnostics
├── module1_group_size.R          Module 1: group size (standalone)
├── module1b_identification.R     Module 1b: homophily vs status maximisation
├── module2_spatial.R             Module 2: spatial confounding
├── module3_mediation.R           Module 3: mediation under confounding
├── run_module4.R                 Module 4: matching mechanisms (crash-safe runner)
├── run_modules1to3.R             Modules 1 and 1b: crash-safe runner
├── results/                      Output CSVs and diagnostic plots
└── module5/                      [Experimental] Boundary estimation comparison
    ├── README.md
    ├── module5_config.R
    ├── module5_dgp.R
    ├── module5_clogit.R
    ├── module5_evaluate.R
    ├── module5_paper2.R
    ├── module5_run.R
    └── results/
```

### utils.R: Shared Infrastructure

All shared functions are in `utils.R`, organised in 17 sections:

| Section | Key functions | Purpose |
|---------|--------------|---------|
| 1 | `generate_agents()` | Create agent populations with education and group membership |
| 2, 3 | `place_agents_spatially()`, `sort_into_workplaces()` | Spatial placement and workplace assignment |
| 4 to 6 | `construct_choice_sets()`, `simulate_choices()` | Choice set construction and partner selection |
| 7, 8 | `estimate_clogit()`, `run_monte_carlo()` | Estimation wrapper and Monte Carlo loop |
| 9 to 13 | `assess_recovery()`, plotting, diagnostics | Bias, RMSE, coverage assessment |
| 14 | `bilateral_matching()`, `calibrate_threshold()` | Bilateral mechanism and threshold calibration |
| 15 | Distance-based encounter functions | For Module 2 spatial matching |
| 16, 17 | `sequential_onesided_matching()`, `deferred_acceptance_matching()` | Alternative matching mechanisms |

## Running the Simulations

### Requirements

- R ≥ 4.0
- Packages: `survival`, `dplyr`, `tidyr`, `ggplot2`, `truncnorm`
- Memory: approximately 4 GB per module at N = 1,000 and R = 200
- Time: 30 minutes to 4 hours per module depending on conditions

### Installation

```r
install.packages(c("survival", "dplyr", "tidyr", "ggplot2", "truncnorm"))
```

### Execution

Each module is self-contained and sources `utils.R`. Modules are independent and can run in any order. All scripts use `set.seed(42)` for reproducibility.

```bash
Rscript run_modules1to3.R module1     # Module 1 (~30 min)
Rscript run_modules1to3.R module1b    # Module 1b (~60 min)
Rscript module2_spatial.R             # Module 2 (~30 min)
Rscript module3_mediation.R           # Module 3 (~3 to 4 hours)
Rscript run_module4.R                 # Module 4 (~2 to 3 hours)
```

For Modules 1 and 1b, the crash-safe runner `run_modules1to3.R` is recommended over the standalone scripts. It uses sequential execution with incremental CSV saves to prevent memory issues (see Technical Notes). Module 4's runner uses the same pattern. Both runners accept a `quick` argument for reduced-scale sanity checks.

### Output

Results are saved to `results/` as CSV files and diagnostic plots:

| File | Contents |
|------|----------|
| `module1_results.csv` | Module 1: oracle + bilateral, 4 compositions |
| `module1b_results.csv` | Module 1b: 3 DGPs × 4 specifications |
| `module2_results.csv` | Module 2: 6 oracle + 3 bilateral conditions |
| `module3_core_stages.csv` | Module 3 Part I: confounding stages A to D, oracle + bilateral |
| `module3_assumption_sweeps.csv` | Module 3 Part II: 8 assumption sweeps (oracle) |
| `module3_interactions.csv` | Module 3 Part III: interaction tests (oracle) |
| `module4_results.csv` | Module 4: all parts including calibrated sweep |
| `module4_pair_diagnostics.csv` | Module 4 Part 5: pair-level utility data |
| `module4_status_diagnostics.csv` | Module 4 Part 6: V distribution diagnostics |
| `module2_*.pdf`, `module3_*.pdf` | Diagnostic plots |

### Reading the CSV Results

Each row in a results CSV represents one parameter under one condition, summarised across Monte Carlo repetitions. The key columns:

| Column | Definition |
|--------|------------|
| `bias` | Mean estimate minus true value |
| `rel_bias_pct` | 100 × bias / \|true value\| |
| `rmse` | Root mean squared error |
| `coverage_95` | Fraction of repetitions where the true value falls in the 95% CI |
| `sd_est` | Standard deviation of estimates across repetitions |
| `mean_pair_rate` | Fraction of agents successfully matched (bilateral/DA only) |

**Note on Module 4 CSV**: `module4_results.csv` has variable column counts because Parts 2 to 4 append extra columns (threshold, acceptance rate, encounter rounds). Read with `fill = TRUE` in R or handle ragged rows in Python.

## Module Details

### Module 1: Group Size and Uneven Trait Distributions

Tests whether clogit recovers homophily and group-preference parameters when education is unevenly distributed across groups of different sizes.

**Conditions**: Four group compositions (balanced 33/33/33, moderate 50/30/20, unbalanced 70/20/10, extreme 85/10/5), each with group-specific education distributions, tested under oracle and bilateral matching.

**Result**: Under oracle, bias is near zero (< 1%) across all conditions. Under bilateral, edu_diff is overestimated by 33 to 47%, roughly constant across compositions.

### Module 1b: Identification of Homophily vs Status Maximisation

Tests whether clogit can distinguish between two preference structures that produce similar matching patterns.

**Design**: Three DGPs (homophily, status_max, mixed) × four clogit specifications. Ego characteristics are conditioned out by clogit strata, so alt_edu and |edu_diff| are not collinear: the absolute value introduces a nonlinearity that permits identification.

**Result**: The correct specification recovers all parameters with < 1.5% bias. Misspecification produces detectable, directional bias. Fitting a homophily model to status_max data yields a small spurious homophily coefficient (+0.127); fitting a status model to homophily data correctly returns near-zero on alt_edu.

### Module 2: Spatial Opportunity vs Homophily

Tests whether residential segregation by education creates spurious homophily estimates when propinquity also drives partner choice.

**Design**: Agents on a 100×100 spatial grid with tunable education-based segregation. Six oracle experiments vary segregation strength, local market structure, and distance measurement resolution. Three bilateral experiments test whether spatial confounding interacts with two-sided matching.

**Result**: Under oracle, the full model (edu_diff + log_distance) recovers both parameters with < 1.5% bias at all segregation levels. Under bilateral, edu_diff is overestimated by about 55% even with distance controlled.

### Module 3: Mediation via Sequential Opportunity Variables

Tests the sequential mediation decomposition used in the dissertation (M0: education only, then M1: + workplace) under unobserved confounding.

**Part I**: Four stages of increasing confounding (no confounder, independent confounder, confounder correlated with education, confounder affecting choice sets), tested under oracle and bilateral matching.

**Part II**: Eight assumption sweeps, varying one structural parameter at a time: pathway strength, choice set size, confounder type, group structure, multiple mediators, workplace concentration, direct effect magnitude, and confounding strength.

**Part III**: Interaction tests (confounding × choice set size, confounding × group structure).

**Core finding**: Under oracle, the education coefficient is robust across all conditions (bias < 5.4%); the workplace coefficient absorbs the confounding; mediation shares are stable but upward-biased. Under bilateral matching, education coefficient bias is about 45% but does not interact with confounding: the biases are additive.

### Module 4: Matching Mechanisms

The central methodological contribution. Tests how the matching mechanism affects clogit parameter recovery.

**Part 1** (Core): Oracle × bilateral × deferred acceptance, under three DGPs.

**Part 2** (Calibrated threshold sweep): Bilateral matching at 5 acceptance rates (10 to 80%) × 2 DGPs. At equivalent market selectivity, homophily is overestimated in selective markets and attenuated in permissive ones, while status maximisation is attenuated at all selectivity levels.

**Part 3** (Group interaction): Bilateral bias across balanced, moderate, and extreme group compositions.

**Part 4** (Rounds sweep): Encounter rounds from 5 to 200, separating encounter limitation from the bilateral bias itself.

**Part 5** (Pair diagnostics): Utility distributions for matched pairs showing the mutual acceptance truncation that drives the bias.

**Part 6** (V distribution diagnostics): Confirms that calibrated thresholds equalize selectivity across DGPs.

## Technical Notes

### macOS Memory and Parallelism

`mclapply` (fork-based parallelism) corrupts R's internal memory state on macOS when used with compiled C code from `survival::clogit`. All production runners use sequential execution (`N_CORES = 1`) with incremental CSV saves after each condition. On Linux, parallelism may work but has not been tested at scale.

### Choice Set Construction

For one-sided matching, each agent's choice set contains J random alternatives from the population. For bilateral matching, after pairs form, each matched agent's choice set contains their actual partner (chosen = 1) plus J minus 1 random alternatives from the full population (chosen = 0). This mirrors the empirical approach where choice sets are sampled from the population regardless of who actually matched.

### Reproducibility

All scripts use `set.seed(42)`. Results may vary slightly across R versions due to differences in random number generation, but qualitative findings should be stable.

## Citation

If using this code, please cite:

> Lindmarker, J. (2026). *Assortative Mating in Sweden: Opportunity, Preference, and the Structure of Partner Choice.* Doctoral dissertation, Linköping University.

## License

MIT License. See LICENSE file for details.
