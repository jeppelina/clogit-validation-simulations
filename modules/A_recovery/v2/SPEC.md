# Module A (version 2): does a well-specified conditional logit do what it aims to do?

**Status:** spec, 2026-10-13, written before the runs. Version 1 is in
`v1/`, and its README says why it was replaced.

## The claim

> A well-specified conditional logit recovers ancestry assortativity free of
> group size, residential segregation and choice-set sampling.
>
> **Well-specified** means:
> - egos from every group are pooled;
> - each ego group gets its own same-group term;
> - each alter group gets a constant;
> - education (and distance, when there is space) is included.
>
> This is the logit form of a log-linear model of the couple table.
>
> - Under one-sided choice it recovers the true taste.
> - Under two-sided matching it recovers the pair-level assortativity, up
>   to a scale set by search frictions (`VALIDATION.md`).

## Setup (unchanged from version 1)

- 1,000 men and 1,000 women; three groups; continuous education.
- Utility: 1 x same group - 1 x abs(edu difference) (- 0.5 x ln distance in
  A2). Identical for both sexes and all groups.
- Settings:
  - `oracle` (one-sided);
  - `gs` (Gale-Shapley, everyone meets everyone);
  - `seq` (sequential search, budget 50).
  Both two-sided settings use 85% pairing.
- Men as egos. Version 1 showed women's estimates match men's in this
  symmetric design.
- 100 uniformly sampled alternatives.
- 200 replications per cell.

## Specifications fitted to every simulated market

| Name | Model | Role |
|---|---|---|
| `full` | same group x ego group (3 terms) + alter-group constants + education (+ distance) | **The test** |
| `no_constants` | same as `full` without alter-group constants | Demonstration: what goes wrong without the column margin |
| `one_group` | majority egos only, alter-group dummies + education (+ distance) | Demonstration: Paper 3's design |
| `loglinear` | Poisson on the 3 x 3 couple table: rows + columns + group diagonals | Reference (A0 and A1 only; it cannot hold distance or sampling) |

## Blocks

| Block | Varies |
|---|---|
| A0 | Baseline: 70/20/10, no space, J = 100 |
| A1 | Group shares 33/33/33, 85/10/5, 90/8/2 |
| A2 | Segregation 0, 0.5, 0.9 (space; `full` with and without distance) |
| A3 | J = 30; stratified sampling (50 of 100 from the ego's own group), with and without the offset (`full` only) |

## Pass criteria for `full`

rho = estimate / target. The target is 1 under oracle and 2 (the pair
surplus) under two-sided matching.

1. **Oracle:** every same-group term and education have abs(rho - 1) < 0.02.
2. **Two-sided:**
   - the three groups' same-group rhos are equal within 0.05;
   - each stays within 0.05 of its A0 value across A1 to A3.
   - The very small group in 90/8/2 (20 men) is reported but not held to
     this, since its estimate is too noisy.
3. **Two-sided, A0 and A1:** `full` agrees with `loglinear` within 0.05.

**Amendment (2026-10-13, after the 3-rep smoke test, before any production
results):**
- Each check passes if the deviation is within its tolerance (0.02 or
  0.05) **or** within 2.5 Monte Carlo standard errors, whichever is larger.
- Reason: the alter-group constants make the small groups' terms noisier.
  A fixed 0.02 bound would fail on simulation noise alone.
- The log-linear agreement check keeps a fixed 0.05: both models are fitted
  to the same markets, so their difference has little noise.

Expected to fail, and reported as demonstrations:
- `no_constants` and `one_group` under two-sided matching;
- `full` without distance in A2;
- stratified sampling without the offset in A3.

## Outputs

`results/v2_A0.csv` to `v2_A3.csv`, `REPORT.md` and one figure. The report leads with
`full`. The demonstrations get one short section.
