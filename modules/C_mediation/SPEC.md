# Module C: mediation under unobserved confounding

**Status:** spec, 2026-10-14, written before the runs. Version 1 is in `v1/`.

## The question

When an unobserved trait drives both who ends up in a setting and whom
they choose, where does the bias land:
- on the **setting** coefficient;
- on the **homogamy** (education) coefficient;
- or on the **mediation share**?

## The toy

- 1,000 men and 1,000 women, one group, education N(0, 1).
  - Groups are not needed here, and the setup is otherwise the same as
    Module A.
- **One kind of setting** (think workplace or university): 10 settings.
  - People join setting k with probability proportional to
    exp(-5 x (abs(edu - c_k) + abs(t - d_k))).
  - c_k is the setting's education level. d_k is its level of the
    unobserved trait.
  - So settings sort on education, which makes them a pathway for
    educational homogamy, and also on t.
- **An unobserved trait t** ("lifestyle"), N(0, 1). It sorts people into
  settings in every cell.
- **Tastes:** -1 x abs(edu difference) + b_z x same setting - gamma x
  abs(t difference).
- **Settings:**
  - one-sided (the benchmark);
  - two-sided (Gale-Shapley, men proposing, 85% pairing).
- **Estimation:** men as egos, 100 uniformly sampled alternatives.

## Cells

- **Confounding strength:** gamma in {0, 0.5, 1}. With gamma = 0, t sorts
  people into settings but does not affect partner choice, so there is
  nothing to confound.
- **Real effect of the setting:** b_z in {0 (placebo), 1 (weaker), 2.5
  (strong)}.
- Both settings (one-sided and two-sided).
- 3 x 3 x 2 = 18 cells, 200 replications each.

## Models

| | Model | Who can fit it |
|---|---|---|
| researcher | M0: education. M1: education + same setting. KHB share from the pair. | A researcher |
| truth | The same, with abs(t difference) added to both | Only the simulation |

**The KHB share.**
- The reduced model uses the setting residualised on education (and on
  abs(t difference) in the truth model).
- share = (b_edu,reduced - b_edu,full) / b_edu,reduced.

## What is reported (per cell, researcher minus truth)

1. **Setting coefficient:** researcher's estimate minus the truth model's.
2. **Education coefficient in M1** (direct homogamy): researcher minus
   truth.
3. **Mediation share:** researcher minus truth.

Also reported against the true tastes under one-sided choice:
- the truth model recovers b_z;
- the truth model recovers the education coefficient, -1.

## Expectations (from version 1 and the archived module 3)

- The confounder lands mostly on the **setting** coefficient: inflated, and
  more so the stronger gamma is.
- The **direct education** coefficient is largely protected.
- The **share** is inflated modestly.
- In the placebo cells, under confounding, a setting that does nothing gets
  a positive coefficient and a positive share.
- Two-sided matching should not change where the bias lands.

No pass/fail criterion: the size and location of the bias are the
findings. Sanity checks:
- with gamma = 0, researcher equals truth within 2.5 Monte Carlo SEs;
- under one-sided choice, the truth model recovers b_z and -1 within
  max(0.05, 2.5 MCSE).

## Outputs

- `results/C.csv`, `summary.csv`
- one figure: bias in each of the three quantities, by gamma, faceted by
  b_z and setting
- `REPORT.md`
