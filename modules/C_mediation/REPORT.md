# Module C report: where unobserved confounding lands in a mediation analysis

**Run:** 2026-10-14. `run_C.R` (200 reps per cell, 18 cells), `analyse_C.R`.
Version 1 is in `v1/`.

**Setup.**
- 1,000 men and 1,000 women; 10 settings (workplaces, say) that people join
  partly by education and partly by an unobserved trait t.
- Tastes: -1 x education difference + b_z x same setting - gamma x
  t difference.
- The researcher fits education only (total homogamy), then education +
  same setting (KHB).
- The truth is the same models with t's difference added, which only the
  simulation can do. Researcher and truth are fitted to the same markets.

## Answer

**Unobserved confounding lands mainly on the setting coefficient.**
- The direct homogamy coefficient moves much less: up to 10% toward zero.
- The mediation share is inflated by a few points. That is large when the
  real channel is weak, modest when it is strong.
- Two-sided matching does not change where the bias lands.

![](results/fig_confounding.pdf)

| Setting effect | gamma | Setting coefficient: researcher / truth (one-sided) | Direct homogamy: researcher / truth | Share: researcher / truth |
|---|---|---|---|---|
| placebo (0) | 0.5 | 0.32 / 0.01 | -0.94 / -1.00 | 0.029 / 0.001 |
| placebo (0) | 1 | **0.56 / 0.00** | -0.90 / -1.00 | **0.052 / 0.000** |
| weaker (1) | 0.5 | 1.32 / 1.00 | -0.95 / -1.00 | 0.111 / 0.083 |
| weaker (1) | 1 | **1.54 / 0.99** | -0.90 / -1.00 | **0.135 / 0.082** |
| strong (2.5) | 0.5 | 2.82 / 2.51 | -0.96 / -1.00 | 0.215 / 0.189 |
| strong (2.5) | 1 | **3.03 / 2.50** | -0.95 / -1.01 | **0.230 / 0.187** |

One-sided choice. Monte Carlo SEs are 0.002 to 0.005.

1. **The setting coefficient absorbs the trait.**
   - The bias is +0.3 at gamma 0.5 and +0.55 at gamma 1, nearly the same
     whatever the setting's real effect.
   - A setting that does nothing gets a coefficient of 0.56.
   - In relative terms, the weaker real effect is overstated by 55% and the
     strong one by 21%.
   - Two-sided matching gives the same bias in units of the homogamy
     coefficient (0.27 and 0.45, against 0.31 and 0.55 one-sided).
2. **The direct homogamy coefficient is only partly protected.** It shrinks
   by 5 to 10% toward zero. Part of the similarity in t that couples share
   is credited to education. Under two-sided matching with a strong setting
   the sign flips slightly: the coefficient is overstated by up to 2%.
3. **The mediation share is inflated by a few points.**
   - +0.03 to +0.05 under one-sided choice, +0.02 to +0.04 two-sided.
   - For a placebo setting that is a share of 0.05 where the truth is 0.
   - For a weaker setting, 0.135 against 0.082, a 64% overstatement.
   - For a strong setting, 0.23 against 0.19, about 20%.

**Checks.**
- The truth model recovers the true setting effect and -1 under one-sided
  choice in every cell (2.49 to 2.51, 0.99 to 1.00, about 0).
- With gamma = 0, researcher and truth agree to within 0.013 on
  coefficients and 0.001 on the share. The SPEC's 2.5-MCSE rule flags a
  few of these, because paired differences have very small MCSEs. The
  sizes are negligible.

## For the papers

- **A setting's coefficient** (workplace, university, neighbourhood) is the
  quantity most exposed to unobserved sorting. Read it as an upper bound on
  the setting's effect when people select into settings on traits that also
  matter for partner choice.
- **The direct homogamy coefficient** is fairly robust: within 10% here.
- **Mediation shares** are inflated by a few points. This matters most when
  the share is small: then a large part of it can be selection rather than
  the setting.
