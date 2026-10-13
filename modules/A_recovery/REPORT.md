# Module A report (version 2): does a well-specified conditional logit do what it aims to do?

**Run:** 2026-10-13. `run_A.R` (200 reps per cell), `analyse_A.R`,
diagnostics `diag_segregation.R`, `explain_group_size.R` and
`explain_loglinear.R`. Version 1 is in `v1/`.

## Answer

**Yes, with one limit.** A conditional logit that pools egos from every
group, gives each ego group its own same-group term, and gives each alter
group a constant:
- recovers the true taste under one-sided choice;
- recovers the same quantity as a log-linear model of the couple table
  under two-sided matching, free of group size and of choice-set sampling.

The limit: **under heavy residential segregation**, two-sided matching
shifts about 10% of the weight from distance to the same-group terms. That
happens even with distance and alter-group constants in the model.

## The well-specified model

    chosen ~ same_group x ego_group (one term per group) + alter-group constants
             + abs(edu difference) [+ ln(distance)] + strata(ego)

**Why the alter-group constants matter.** In a two-sided market each group's
members differ in how attainable they are as partners. Each mixed couple is
an out-marriage for one person from each group, so mixing is shared. The
constants absorb that, as the column margin does in a log-linear model.
Without them, the same-group terms absorb it and move with group size
(figure, bottom row; version 1).

## Results

rho = estimate / target. The target is the true taste (1) under one-sided
choice, and the pair surplus (2) under two-sided matching. Settings:
- one-sided (oracle);
- Gale-Shapley (everyone meets everyone);
- sequential search (about 17 meetings).

Both two-sided settings use 85% pairing. Monte Carlo SEs are 0.01 to 0.02,
and up to 0.04 for the 5% group.

![](results/v2_figure.pdf)

**1. Group size (A0, A1).** Same-group rho by group: largest / middle /
smallest.

| Shares | One-sided | Gale-Shapley | Sequential search |
|---|---|---|---|
| 33/33/33 | 1.01 / 1.00 / 1.01 | 0.96 / 0.96 / 0.96 | 0.65 / 0.65 / 0.65 |
| 70/20/10 | 1.01 / 1.00 / 1.02 | 0.94 / 0.97 / 0.98 | 0.66 / 0.65 / 0.62 |
| 85/10/5 | 0.93 / 1.07 / 0.95 | 0.97 / 0.97 / 0.95 | 0.64 / 0.69 / 0.66 |

- **One-sided:** the true taste, at every composition.
- **Gale-Shapley:** the pair surplus at a scale of about 0.95, equal across
  groups and compositions.
- **Sequential search:** a scale of about 0.65, nearly equal across groups.
  The smallest group sits up to 0.04 lower, and the log-linear model shows
  the same (0.663 / 0.650 / 0.621 at 70/20/10). So this is a feature of the
  frictional market, not of the logit: with about 17 meetings, a member of
  a 10% group meets fewer than two people from their own group.
- **The logit agrees with the log-linear model within 0.01** in every cell
  of A0 and A1 (as checked, apart from 90/8/2).
- **90/8/2 is underpowered and left out.** The smallest group has about 15
  egos. In 27% of markets its term separates, because every member pairs
  within the group. Medians are close to target (Gale-Shapley 0.93 / 1.03 /
  0.91), means are not. The log-linear model breaks down in the same way
  (empty cells).

**2. Sampling (A3).**

| Design | One-sided | Gale-Shapley | Sequential search |
|---|---|---|---|
| J = 100, uniform | 1.02 / 0.97 / 0.99 | 0.95 / 0.97 / 0.97 | 0.68 / 0.65 / 0.62 |
| J = 30, uniform | 1.02 / 0.97 / 0.98 | 0.95 / 0.97 / 0.97 | 0.68 / 0.65 / 0.62 |
| Stratified (50 of 100 from own group), with offset | 1.02 / 0.97 / 0.99 | 0.95 / 0.96 / 0.97 | 0.68 / 0.65 / 0.62 |
| Stratified, **without** offset | 1.87 / -0.42 / -1.21 | 1.38 / 0.27 / -0.13 | 1.10 / -0.04 / -0.48 |

- Fewer alternatives and stratified sampling change nothing, **as long as
  the McFadden-Manski offset is used**.
- Without the offset, the alter-group constants absorb the oversampling and
  the same-group terms are badly wrong.

**3. Segregation (A2).** Distance is in the model.

| Segregation | Setting | Same-group rho (g1 / g2 / g3) | Distance rho |
|---|---|---|---|
| 0 | oracle | 0.97 / 1.04 / 1.04 | 1.01 |
| 0.9 | oracle | 1.00 / 0.99 / 0.99 | 1.00 |
| 0 | Gale-Shapley | 0.95 / 0.95 / 0.94 | 0.93 |
| 0.9 | Gale-Shapley | **0.99 / 1.05 / 1.05** | **0.81** |
| 0 | sequential | 0.66 / 0.63 / 0.61 | 0.57 |
| 0.9 | sequential | **0.76 / 0.74 / 0.71** | **0.47** |

- Under one-sided choice, distance fully absorbs spatial endogamy.
  Without distance, the same-group terms rise to 1.56 to 1.66.
- **Under two-sided matching with heavy segregation, the distance
  coefficient loses about 0.1 of its scale and the same-group terms gain
  0.05 to 0.10.**
- Moderate segregation (0.5) shows no such shift.
- Two explanations are ruled out:
  - a mis-specified distance gradient (a 5-df spline did not help, v1);
  - local group margins (alter-group constants by region did not help,
    `diag_segregation.R`).
- The likely reason: two-sided matching does not produce an exactly logit
  pattern. When two regressors are strongly correlated, as group and
  distance are under heavy segregation, the misfit is not shared evenly
  between them.
- Consequence: with very strong segregation, residual endogamy net of
  distance is modestly overstated, by about 10%. The education term is
  unaffected (0.95 to 0.96).

## Demonstrations (expected to fail)

- **Without alter-group constants:** same-group terms move with group size.
  Under Gale-Shapley at 85/10/5: 0.56 for the largest group, 1.30 and 1.36
  for the others (figure, bottom row).
- **One ego group only (Paper 3's design).** The majority's terms against
  the minority groups fall as the majority grows. Under Gale-Shapley:
  - 33/33/33: 0.97 / 0.97;
  - 70/20/10: 0.73 / 0.64;
  - 85/10/5: 0.59 / 0.52.

  The constants and the same-group term cannot be separated with one row
  of the couple table. This is the source of Module D's spurious
  hardening.
- **No distance under segregation; stratified sampling without the offset:**
  both fail as expected (tables above).

## What this means for the papers

- **Fit pooled models with alter-group constants**, or equivalently a
  log-linear model of the couple table, to compare assortativity across
  groups of different size, or over time when composition changes.
- **Use the offset whenever choice sets are stratified.**
- **Read coefficient levels as pair-level assortativity scaled by market
  frictions** (about 0.95 to 0.65 here), not as individual tastes. Use them
  to compare groups and terms within the same market.
- **Under heavy segregation**, treat residual endogamy net of distance as
  modestly overstated.
