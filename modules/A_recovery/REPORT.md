# Module A report (version 3): four minimal tests

**Run:** 2026-10-14. `run_A.R` (200 reps per cell), `analyse_A.R`,
`plot_layouts.R`.

**Setup.**
- Two groups.
- One-sided choice (what the logit assumes) vs Gale-Shapley (both must
  agree, 85% pair).
- Well-specified model: same group + minority-alter constant + education.

**Targets.**
- One-sided: the true taste.
- Two-sided: the pair surplus (the sum of both sides' terms). It is
  recovered at a scale of about 0.96 when everyone meets everyone.

**All 28 pre-registered checks pass** (`results/checks.csv`). Monte Carlo
SEs are 0.005 to 0.02.

## Test 1. Group size: the alter-group constant makes the estimate free of group size

![](results/fig_test1_group_size.pdf)

| Minority share | One-sided, with constant | Two-sided, with constant | Two-sided, without constant |
|---|---|---|---|
| 50% | 1.00 | 0.97 | 0.97 |
| 20% | 1.01 | 0.97 | 0.98 |
| 5% | 0.96 | 0.96 | **0.82** |

Same-group estimate divided by its target.

- With the constant, the estimate does not move with group size under
  either setting.
- Without it, a two-sided market makes a small minority look less
  assortative than it is (0.82). The same couples are still there; the
  constant soaks up how attainable minority partners are.

## Test 2. Segregation: with ln(distance), both layouts are recovered

![](results/layouts.pdf)
![](results/fig_test2_segregation.pdf)

| Layout | One-sided, with distance | One-sided, without | Two-sided, with distance | Two-sided, without |
|---|---|---|---|---|
| mixed | 1.00 | 1.00 | 0.96 | 0.95 |
| two sides, moderate | 1.01 | 1.09 | 0.96 | 1.03 |
| two sides, strong | 1.00 | **1.25** | 0.96 | **1.21** |
| enclaves, moderate | 1.00 | 1.04 | 0.96 | 1.00 |
| enclaves, strong | 1.02 | **1.14** | 0.96 | **1.11** |

Same-group estimate divided by its target. The distance term is recovered
too: 0.98 to 1.01 one-sided, 0.92 two-sided.

- **With distance in the model, segregation does not affect the
  same-group estimate in either layout or setting.**
- **Without it, propinquity is read as endogamy**: up to 25% too high when
  each group has its own side, up to 14% with enclaves.
- The two-sides layout inflates more at similar levels of segregation:
  minority exposure 0.70 against 0.78, dissimilarity 0.91 against 0.85.
- Module A v2 found that, under heavy segregation, two-sided matching
  shifted about 10% of the weight from distance to the same-group terms.
  **That does not recur here.** The v2 layout packed each of three groups
  tightly around a single point (sd 0.08), a far more extreme
  concentration than either layout here. The shift therefore seems to need
  extreme clustering. The v2 report keeps the record.

## Test 3. Sex differences: two-sided data cannot say which sex holds a taste

![](results/fig_test3_sex.pdf)

| World (men's, women's same-group taste) | One-sided, men / women | Two-sided, men / women |
|---|---|---|
| only men care (2, 0) | 1.99 / -0.01 | 1.85 / 1.85 |
| only women care (0, 2) | 0.00 / 2.02 | 1.85 / 1.85 |
| both care (1, 1) | 1.00 / 0.99 | 1.91 / 1.91 |
| minority 60% male (1, 1) | 1.00 / 1.00 | 1.94 / 1.94 |

Estimates, not ratios.

- One-sided choice recovers each sex's own taste.
- **Two-sided: the men's and women's models give the same number in every
  world** (differences of 0.001 or less). They estimate the pair, so they
  cannot say who holds the taste.
- The level depends a little on how the taste is split: 1.85 when one side
  holds it all, 1.91 when it is shared.
- **An unbalanced sex ratio no longer produces a sex gap once the
  minority-alter constant is in the model** (1.94 / 1.94). Module B found
  1.97 / 1.58 without it. That gap was the same margin problem as in
  test 1.

## Test 4. Status vs homophily: status is nearly invisible, and what it leaves looks like homophily

![](results/fig_test4_status.pdf)

| True tastes | One-sided: homophily term / status term | Two-sided: homophily term / status term | Couples' education correlation, two-sided |
|---|---|---|---|
| homophily (-1, 0) | -1.00 / 0.00 | -1.93 / 0.00 | 0.78 |
| status (0, +2) | -0.01 / 2.01 | **-0.56 / 0.22** | 0.35 |
| both (-1, +2) | -1.00 / 2.01 | -2.07 / **0.24** | 0.78 |

Targets: homophily -1 one-sided and -2 two-sided; status 2 in both.

- One-sided choice separates the two.
- **Two-sided: a status taste of 2 shows up as 0.22**, about a tenth.
- When everyone wants the more educated, competition sorts couples by
  education (correlation 0.35). The model reports that sorting as
  homophily: -0.56 in a world with no homophily.
- With both tastes present, homophily is recovered in full and status is
  nearly lost. A researcher would conclude "homophily, little status".
- Terms that describe the alter alone, such as the partner's education
  here and the alter-group constant in test 1, measure who is attainable,
  not what the chooser wants.

## In one paragraph

A well-specified conditional logit does what it aims to do on group size,
segregation and choice-set sampling, as long as distance is in the model
and alter-group constants absorb the margins. Under one-sided choice it
recovers tastes. Under two-sided matching it recovers the pair-level
assortativity, at a stable scale of about 0.96. What it cannot do under
two-sided matching is attribute that assortativity:
- to one sex rather than the other (test 3);
- to a taste for similarity rather than competition for status (test 4).
