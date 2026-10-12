# Module D: Decomposing change over cohorts

**Status:** spec, 2026-10-13. Uses the market and settings validated in A.

## The claim under test

> When the partner market changes across cohorts, the composition/coefficient
> split (Paper 3) attributes the change to the right component. In
> particular:
> - composition change alone produces no coefficient component;
> - a real change in tastes produces a coefficient component of the right
>   sign.

Module A gives reason to doubt the first part under two-sided matching. The
same-group coefficient loses scale as the majority grows (A1). A market
in which the majority shrinks would then show a coefficient change with no
change in tastes. D tests this directly.

## Setup

- Three cohorts, each its own market of 1,000 men and 1,000 women.
- The minorities grow across cohorts. Group shares: 90/7/3, then 80/13/7,
  then 70/20/10.
- Utility: b_g(cohort) x same group - 1 x abs(edu difference), for both
  sexes and all groups.
- Egos: majority-group men. One conditional logit per cohort, with 100
  alternatives drawn uniformly. Separate models per cohort are equivalent to
  Paper 3's full cohort interaction.
- Outcome: the share of majority men's unions that are with majority women.

**Decomposition** (cohort 1 to cohort 3, as in Paper 3). P(data, b) is the
predicted majority share: the mean over egos of the predicted probability
on same-group alternatives.

- composition = P(data_3, b_1) - P(data_1, b_1)
- coefficient = P(data_3, b_3) - P(data_3, b_1)
- composition + coefficient = the observed change. The clogit score
  equation makes P(data_c, b_c) equal the observed share exactly.

**The truth, by counterfactual.**

- Cohort 3's population is matched a second time with cohort 1's tastes,
  using the same people, shocks and meetings.
- true coefficient = endogamy(cohort 3, own tastes) - endogamy(cohort 3,
  cohort 1's tastes)
- true composition = endogamy(cohort 3, cohort 1's tastes) - endogamy(cohort
  1, own tastes)
- This works under any matching rule, including two-sided ones where
  per-side tastes are not identified.

## Cells

| Scenario | b_g by cohort | Settings |
|---|---|---|
| null | 1, 1, 1 | oracle, gs, seq |
| harden | 1, 1.25, 1.5 | oracle, gs, seq |
| soften | 1, 0.75, 0.5 | oracle, gs, seq |
| frictions | 1, 1, 1; search budget 25, 50, 100 by cohort | seq only |

- The outside option is calibrated to 85% pairing for every cohort and
  counterfactual market.
- 200 replications per cell.

## Pass / fail

- **null:** the estimated coefficient component is within 2 MC standard
  errors of zero, and below 0.002 in absolute value, in every setting.
- **harden / soften:** the estimated coefficient component has the sign of
  the true component in at least 95% of reps, and its mean is within 20% of
  the mean true component.
- **frictions:** reported. The true taste change is zero, so any coefficient
  component shows that the component measures more than tastes.
