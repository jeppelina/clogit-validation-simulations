# Module A: Recovery

**Status:** draft spec for sign-off, 2026-10-12. No code written yet.
**Replaces:** old Module 1 (group size), Module 2 (spatial), the oracle half of
Module 1b, and the selectivity sweep in Module 4.

## The one claim

> Under the conditions the papers actually face (unequal group sizes,
> residential segregation, county choice sets, sampled and stratified
> alternatives), conditional logit recovers its target. Under one-sided
> (oracle) choice the target is the chooser's own coefficient. Under stable
> two-sided matching (Gale-Shapley) the target is the pair surplus, the sum of
> both sides' coefficients.

A has a second job: it is the control condition for B, C and D. If A
passes, any distortion found later comes from the departure that module
introduces, not from the pipeline.

## Why the Gale-Shapley target is the pair surplus

When both partners must prefer the match to their alternatives, the
frequency of a pair type depends on what both sides get from it. In large
markets with logit taste shocks, stable matching produces match frequencies
in which the two sides' systematic utilities enter as a sum (Dagsvik 2000;
Menzel 2015), which is the Choo-Siow (2006) pair surplus. So a one-sided
clogit fitted to stable matches should estimate beta_M + beta_W, not
beta_M. The April appendix measured two-sided bias against beta_M alone,
which is part of why it reported "overestimation".

Whether this limit holds at our market sizes is not assumed. Block A0 tests
it. Two acceptable outcomes:

- **Recovers:** estimate / surplus is within 5% of 1.
- **Stable bias:** the ratio departs from 1 but holds steady across the
  A1 to A3 conditions. We then report the ratio as a known scale factor.

The failure that matters is a ratio that **moves with the conditions**: that
would mean group size, segregation or sampling interacts with
two-sidedness, and the papers' comparisons across groups or designs would
inherit it.

## Population and utility (shared by every block)

- Two sexes, N per sex (default 1,000; A0 varies it).
- Three ancestry groups (shares vary in A1; default 70/20/10).
- Education in three levels, to mirror the papers' `Edudiff3`.
- Location on a 2D grid divided into 9 counties; groups concentrated
  according to a segregation parameter (A2).
- Utility of man m for woman w, and of w for m:

      U_mw = b_g * same_group + b_e * same_edu + b_d * ln(dist) + eps_mw
      V_wm = b_g * same_group + b_e * same_edu + b_d * ln(dist) + eta_wm

  All three terms are symmetric in the pair. **Both sexes have identical
  coefficients in A**; asymmetry between the sexes is Module B's subject.
  So the oracle target is b, and the Gale-Shapley target is 2b.
- Shocks eps and eta are i.i.d. Gumbel, **drawn once per pair and held
  fixed** for the whole matching process.
- Each agent has an outside option, u0 + Gumbel shock. Partners below it are
  unacceptable. u0 is calibrated to a target pairing rate.
- Default coefficients (per side): b_g = 1.0, b_e = 0.5, b_d = -0.5. To be
  confirmed by a calibration gate (see below), not taken as given.

## Matching

- **Oracle:** each ego picks the best option, including staying single, from
  the full opposite-sex pool. No capacity limit. This is exactly the process
  clogit assumes.
- **Gale-Shapley, men proposing; Gale-Shapley, women proposing.**
  - Full preference lists, truncated at the outside option.
  - Run to convergence, with no round cap.
  - Verified stable: no blocking pairs, checked exhaustively on small
    markets in `tests/`.
- Threshold search is not used.

**Neither existing implementation qualifies.**

- `utils.R::deferred_acceptance_matching` runs on a single-sex pool.
- Module 7's deferred acceptance, also called as `gale_shapley` by Module 8,
  has three problems:
  - each proposer samples 8 random candidates per round;
  - rounds are capped;
  - shocks are redrawn at every comparison, so preferences are not fixed and
    the outcome is not a stable matching.

A uses a new implementation in `R/matching.R`. Separately, the audit has to
record that Module 7 Exp 6 and Module 8's "Gale-Shapley" results come from
that limited-search process (see "Consequences outside A").

## Estimation

Mirrors the papers:

- Separate models with men as egos and women as egos.
- Egos are partnered individuals only.
- The actual partner plus J - 1 alternatives sampled from the opposite-sex
  pool.
- Estimator: `survival::clogit`.

By symmetry the men's and women's estimates should agree in A. That
agreement is a built-in check on the code.

## Blocks

| Block | Question | Varies | Matching |
|---|---|---|---|
| **A0 Anchor** | With no threats present, what does each mechanism recover? | N per sex {500, 1,000, 2,000} x pairing rate {~95%, ~75%, ~50%} | Oracle, GS-men, GS-women |
| **A1 Group size** | Does recovery hold when one group dominates and one is tiny? | Shares {33/33/33, 70/20/10, 85/10/5, 90/8/2} | Oracle, GS |
| **A2 Space** | Does log-distance absorb spatial confounding, and do county choice sets change that? | Segregation {none, moderate, high} x model {with, without ln(dist)} x choice set {national, county} | Oracle, GS |
| **A3 Sampling** | Do the papers' choice-set designs recover the target? | J {30, 100}; uniform vs stratified (oversampled same-group alternatives) with and without the McFadden-Manski offset | Oracle, GS |

**Two rules for the GS column in A1 to A3:**

- GS-women runs only if A0 shows that the proposing side changes the
  estimates.
- Pairing rate is held at the A0 value closest to the papers' setting.

**A3 also carries a kernel check.** On the stratified-with-offset data, fit
the same model with the patched fastclogit kernel (as used in Paper 4) and
with `survival::clogit`. Coefficients should agree to 1e-6, and the
fastclogit fit should converge on the gradient criterion, not on the
plateau rule.

**Calibration gate (before any block runs).** With the default
coefficients, one realisation per mechanism must show:

- majority-group endogamy in a plausible range (80 to 95%);
- no complete separation in any per-sex model;
- the pairing rate at its target.

If not, adjust the coefficients and record the change here before running
the blocks.

**Replications:** 200 per cell. The runtime estimate comes after the GS
implementation is benchmarked.

## Metrics and pass criteria

For each coefficient, cell and ego sex, report:

- mean estimate, relative bias, RMSE, and 95% coverage against the
  mechanism's target;
- under GS, the ratio rho = estimate / surplus.

| Mechanism | Passes if | Fails if |
|---|---|---|
| Oracle | abs(relative bias) < 2% and coverage in [0.93, 0.97] in every cell | Bias appears with small groups, county sets or stratified sampling |
| GS | rho within 5% of 1, **or** rho stable: within 5 percentage points of its A0 value at the same N and pairing rate, across all A1 to A3 cells | rho moves with group shares, segregation or sampling design |

**What the omitted-distance model is for.** In A2, the model without
ln(dist) is expected to fail under segregation. That is a demonstration,
not a test of the pipeline. It is reported as the size of the omitted-
variable bias, under each mechanism.

## Outputs

- `results/A0.csv` to `results/A3.csv`: one row per cell x coefficient x
  ego sex.
- `REPORT.md`: the claim, a pass/fail line per block, and one figure
  (rho across all cells, grouped by block).
- **A comparison table against the old Modules 1, 2 and 4**, so every
  changed number has a stated reason.
- Every results CSV carries the git commit and seed it was produced with.

## Not in A

- Sexes with different coefficients, and which side holds the assortativity
  (Module B).
- Homophily versus status (B).
- Mediation and channel shares (C).
- Change over time (D).
- Threshold search and selectivity sweeps (dropped).

## Consequences outside A

The audit has to answer whether Module 8's Problem B result, and Module 7
Exp 6, survive a real stable matching. In particular: "the sex gap appears
only under competitive matching". Those results currently rest on a
limited-search process with redrawn shocks. Module B will rerun the
relevant cases on the new implementation.

## Open questions for Jesper

1. **County choice sets.** My default:
   - egos are kept regardless of where their partner lives;
   - the actual partner is always in the set;
   - alternatives are drawn from the ego's county.

   Is that what Papers 1 and 3 do with partners from outside the county?
2. **Default J.** I propose 100 (Papers 1 and 3), with 30 as the A3
   variation. Or the reverse?
3. **Education as a three-level same/different term** to mirror `Edudiff3`,
   rather than the old continuous |edu difference|. Agree?
4. **Calibration targets.** Is 80 to 95% majority endogamy and a 75%
   pairing rate a fair stylised Sweden, or do you want them pinned to
   numbers from the papers?
