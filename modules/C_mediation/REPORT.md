# Module C report: mediation shares

**Run:** 2026-10-13. `run_C.R` and `run_C_coarse.R` (200 reps per cell, men
as egos), `analyse_C.R`. KHB share = (b_reduced - b_full) / b_reduced, with
the mediator residualised on focal and controls. Truth = the oracle share
in the same world; for C3, the share with the unobserved trait controlled.

## Summary

| Cell | oracle | gs | seq |
|---|---|---|---|
| C1 placebo channel (true effect 0), continuous education | **-0.001** | 0.004 | 0.005 |
| C2 real channel, continuous education | 0.224 | **0.275** | 0.208 |
| C3 channel with an unobserved confounder (truth in brackets) | 0.231 (0.188) | 0.242 (0.219) | 0.251 (0.208) |
| C4 distance channel, segregation 0.5 (Paper 1) | 0.150 | 0.142 | 0.129 |
| C4 distance channel, segregation 0.9 (Paper 1) | **0.368** | **0.314** | **0.270** |
| Placebo channel, coarse education (same tercile) | **0.083** | 0.084 | 0.082 |
| Real channel, coarse education | **0.524** | **0.714** | 0.499 |

Monte Carlo SEs are 0.0004 to 0.004.

## What holds and what does not

**1. KHB is clean when the homogamy term is right.**
- With continuous education in the model, a channel that does nothing gets
  a share of zero under one-sided choice, and at most 0.005 under two-sided
  matching.
- There is no general "artifact floor" from conditioning on a mediator that
  is correlated with the focal variable. The channel here correlates -0.31
  with abs(edu difference), and the share is still zero.

**2. A coarse homogamy measure creates one.**
- The tastes are continuous in education. When the model measures
  homogamy as "same tercile", an education-sorted channel with no effect of
  its own is credited with 8% of the homogamy coefficient.
- A real channel's share more than doubles: 0.22 becomes 0.52 under
  oracle, and 0.71 under Gale-Shapley.
- The channel absorbs the similarity the coarse measure leaves out. This is
  the mechanism behind the earlier module 8 "floor" (archived), which was
  measured on a coarse "both in the top quarter" indicator. **For Paper 4,
  whose homogamy term is three-level (`Edudiff3`), the risk is real.** The
  fix is a flexible homogamy specification, not subtracting a floor.

**3. Two-sided matching moves the shares, in a direction set by the
market.**
- Institution channel: Gale-Shapley overstates the share (0.275 against
  0.224). Sequential search understates it slightly (0.208).
- The share is a ratio and would survive a common rescaling. It moves
  because two-sided matching rescales the focal coefficient and the
  mediator by different amounts: under gs the focal term keeps about 0.77
  of its pair surplus, the mediator about 0.95.

**4. Two-sided matching understates propinquity's share of endogamy
(Paper 1).**
- At segregation 0.9 the share falls from 0.37 (oracle) to 0.31 (gs) and
  0.27 (seq): 15% and 27% less.
- This is Module A's finding seen from the other side. Under two-sided
  matching, distance does not absorb all of the spatial endogamy, so less
  of it is credited to distance and more stays in the residual endogamy.
- At moderate segregation the gap is small (0.15 against 0.13 to 0.14).

**5. An unobserved trait that sorts people into channels inflates the
share, but modestly here.**
- Inflation is +0.02 to +0.04, about 10 to 20% of the share.
- The size depends on how strongly the trait drives both channel membership
  and partner choice. Here both were set equal to education's role, so this
  is an illustration, not a bound.

## For the papers

- **Paper 4.** The mediated shares are most at risk from the coarse
  homogamy measure (point 2), more than from two-sided matching or
  confounding. Checking the decomposition with a finer education
  specification is the single most informative robustness check.
- **Paper 1.** The share of endogamy attributed to propinquity is, if
  anything, understated under two-sided matching. The paper's conclusion
  that propinquity explains a minority of endogamy is conservative in that
  respect.
