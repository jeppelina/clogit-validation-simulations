# Module D: change over cohorts

**Status:** spec, 2026-10-14, written before the runs. Version 1 is in `v1/`.
Jesper does not expect to use this module, but asked for it to be done.

## The question

When the minority grows between two cohorts, does the well-specified model
separate a change in assortativity from the change in composition?

## The toy

- **Two cohorts.** Each is its own market with the Module A setup: 1,000
  men and 1,000 women, two groups, education N(0, 1).
- **The minority grows from 10% to 30%.**
- **Tastes:** b_g(cohort) x same group - 1 x abs(edu difference), for
  everyone.
- **Three scenarios:**

  | Scenario | b_g, cohort 1 | b_g, cohort 2 |
  |---|---|---|
  | constant | 1 | 1 |
  | harden | 1 | 1.5 |
  | soften | 1 | 0.5 |

- **Settings:**
  - one-sided;
  - two-sided (Gale-Shapley, 85% pairing in every market).
- 3 x 2 = 6 cells, 200 replications each.

## Models

Fitted per cohort, men as egos, both groups pooled, 100 uniformly sampled
alternatives:
- **well-specified:** same group + minority-alter constant + education (as
  in Module A);
- **without constant** (demonstration): same group + education.

## The decomposition

Outcome: the endogamy share, i.e. the share of men's unions within their
own group (one minus the mixed share).

- P(d, b) is the predicted endogamy share using cohort d's choice sets
  and coefficients b. At the MLE, P(d_c, b_c) equals the observed share.
- **composition** = P(d_2, b_1) - P(d_1, b_1)
- **coefficients** = P(d_2, b_2) - P(d_2, b_1). It is split by term
  (same group, minority-alter constant, education) using Shapley values,
  i.e. averaged over every order in which the terms are switched from
  cohort-1 to cohort-2 values.
- The **same-group part** is the decomposition's estimate of the change due
  to assortativity.

## The truth

- Cohort 2's people are rematched with cohort 1's tastes, using the same
  random numbers (and an outside option recalibrated to 85% pairing).
- **true assortativity effect** = endogamy(cohort 2, own tastes) -
  endogamy(cohort 2, cohort 1's tastes). It is exactly 0 in the constant
  scenario.

## Expected

**Well-specified:**
- the same-group part tracks the true effect: zero in the constant
  scenario, the right sign and roughly the right size otherwise;
- the constant's part absorbs the market's response to composition;
- the change in the same-group coefficient itself is zero in the constant
  scenario.

**Without the constant:** the same-group part picks up spurious change
under two-sided matching (Module A test 1).

## Checks

1. **Well-specified, constant scenario, both settings:** the same-group
   part and the change in the same-group coefficient are within max(0.005,
   2.5 MCSE) of 0.
2. **Well-specified, harden and soften:** the same-group part has the sign
   of the true effect in at least 95% of reps. Its mean and the true mean
   are reported side by side.

## Outputs

- `results/D.csv`, `summary.csv`
- one figure: estimated vs true assortativity effect by scenario, setting
  and model
- `REPORT.md`
