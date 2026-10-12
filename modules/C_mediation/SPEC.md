# Module C: Mediation

**Status:** spec, 2026-10-13. Uses the market and settings validated in A.

## The claim

> Mediation shares from KHB decompositions of a conditional logit (Papers 1
> and 4) are ratios of coefficients, so a common rescaling would leave them
> untouched. Two things still distort them:
> - an unobserved trait that drives both channel membership and partner
>   choice inflates the share;
> - two-sided matching rescales the focal coefficient and the mediator by
>   different amounts, which moves the share.
>
> A channel that does nothing should get a share of zero. We check that it
> does.

## Setup

Same market as A (1,000 per sex, 70/20/10, continuous education).

**Settings:** `oracle` (the benchmark: what a one-sided world shows), `gs`
and `seq` (see `R/run_helpers.R`).

**KHB.** Focal variable X, mediator Z, controls C.

- Full model: X + Z + C.
- Reduced model: X + Z_res + C, where Z_res is Z residualised on X and C
  by OLS across all choice-set rows.
- Share = (b_X,reduced - b_X,full) / b_X,reduced.

Residualising keeps the two models on the same scale, which is the point of
KHB.

**Truth.** The oracle share in the same world; for C3, the share with the
trait controlled (infeasible for a researcher, available in a simulation). Tastes are identical across
settings, so a one-sided market's share is what the decomposition should
recover.

## Cells

**Institution channel (Paper 4)**
- Focal: abs(edu difference). Mediator: same institution. Control: same group.
- 10 institutions, each with a characteristic education level. People join
  institution k with probability proportional to
  exp(-s * abs(edu - c_k)), with sorting strength s = 5. That gives about
  10% same-institution overlap among alternatives, and a correlation of
  -0.3 between same institution and abs(edu difference).
- Utility: -1 x abs(edu difference) + 1 x same group + b_z x same
  institution, for both sexes.
- b_z = 2.5 gives an oracle share of about 0.24. These values were chosen so
  the channel is visible, not to match any real data.

| Cell | b_z | Unobserved trait | Question |
|---|---|---|---|
| C1 placebo | 0 | none | Does a channel that does nothing get a zero share? |
| C2 channel | 2.5 | none | Does two-sided matching move the share? |
| C3 confounded | 2.5 | t ~ N(0,1): it sorts people into institutions alongside education, and enters utility as -1 x abs(t difference) | How much does an unobserved trait inflate the share? |

**Distance channel (Paper 1)**
- Focal: same group. Mediator: ln(distance). Control: abs(edu difference).
- Population and utility as in A2 (b_d = -0.5).

| Cell | Segregation | Question |
|---|---|---|
| C4 moderate | 0.5 | Does two-sided matching change the share of endogamy attributed to propinquity? |
| C4 high | 0.9 | Same, under heavy segregation |

A2 already hints at the answer. The naive (non-KHB) share is 0.37 under
oracle, 0.32 under `gs` and 0.28 under `seq` at segregation 0.9.

**Replications:** 200 per cell, men as egos. In A and B, women's estimates
matched men's in every symmetric world.

## Pass / fail

- **C1:** the share is within 0.02 of zero in every setting.
- **C2 and C4:** the two-sided share is within 0.03 of the oracle share.
  Otherwise we report the direction and size of the gap.
- **C3:** we report the inflation under each setting. No pass criterion:
  the inflation is the finding.
