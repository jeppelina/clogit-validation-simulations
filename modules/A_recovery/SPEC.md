# Module A (version 3): four minimal tests of a well-specified conditional logit

**Status:** spec, 2026-10-14, written before the runs. Earlier versions are
in `v1/` and `v2/`. This version also absorbs Module B (sex differences;
status vs homophily).

## The toy

- 1,000 men and 1,000 women. Two groups: a majority and a minority (20%
  unless a test varies it). One continuous education score, N(0, 1).
- Tastes, identical for everyone unless a test varies them:
  +1 x same group - 1 x abs(edu difference).
- **Two settings:**
  - **one-sided:** each person picks freely. This is what the logit
    assumes; it is the benchmark.
  - **two-sided:** Gale-Shapley, men proposing. Everyone meets everyone;
    both partners must prefer the union to staying single. The outside
    option is calibrated so 85% pair up.
- **Estimation:**
  - 100 alternatives sampled uniformly from the other sex;
  - one conditional logit, men as egos (both sexes in test 3);
  - egos from both groups pooled.
- **The well-specified model:**

      same group + minority-alter constant + abs(edu difference)

  plus one term per test.
- With two groups the couple table has one degree of association, so there
  is one same-group coefficient (half the couple-table log odds ratio). The
  minority-alter constant is the column margin: how attainable minority
  partners are.

**Targets.**
- One-sided: the true taste.
- Two-sided: the pair surplus, the sum of both sides' terms (2 for same
  group and education). It is recovered up to a scale near 0.95 when
  everyone meets everyone (Module A v2).
- rho = estimate / target.

**Replications:** 200 per cell.

## Test 1: Group size

- **World:** minority share 50%, 20%, 5%.
- **Models:** well-specified vs **without** the minority-alter constant.
- **Claim:** with the constant, same-group rho does not depend on group
  size, under either setting. Without it, two-sided rho drifts.
- **Pass:** the well-specified rho at 20% and 5% is within max(0.05,
  2.5 MCSE) of its value at 50%. One-sided rho is within that bound of 1.

## Test 2: Segregation

- **World:** propinquity added to tastes (-0.5 x ln distance, both sides).
  Five layouts, plotted in `results/layouts.pdf`:
  - mixed;
  - two sides, moderate and strong: each group mostly on its own side, the
    side's width proportional to the group's share;
  - enclaves, moderate and strong: the minority mostly in three small
    clusters, the majority mostly avoiding them.

  The moderate and strong levels reach similar exposure and dissimilarity
  in both layouts (`results/layouts_stats.csv`).
- **Models:** well-specified + ln(distance), vs well-specified without it.
- **Claim:**
  - with distance, same-group and distance are recovered (one-sided) or
    close to it (two-sided);
  - without distance, the same-group term absorbs propinquity, more so the
    stronger the segregation.
- **Pass (with distance):**
  - one-sided: same-group rho within the bound of 1;
  - two-sided: within max(0.05, 2.5 MCSE) of the mixed layout.
  - Module A v2 found a shift of about 0.1 under heavy segregation. If it
    recurs, it is reported as a limit.

## Test 3: Sex differences

- **Worlds:**

  | World | Men's same-group taste | Women's same-group taste |
  |---|---|---|
  | only men care | 2 | 0 |
  | only women care | 0 | 2 |
  | both care | 1 | 1 |
  | sex ratio | 1 | 1 |

  - The first three worlds have the same pair surplus.
  - In the sex-ratio world the minority is 60% male: 24% of men and 16% of
    women are minority.
- **Model:** well-specified, fitted separately for men and for women.
- **Claim:**
  - one-sided recovers each sex's own taste;
  - two-sided gives the same estimate for men and women in the first three
    worlds, so the data cannot say who holds the taste.
  - The sex-ratio world tests whether the minority-alter constant also
    removes the sex gap that version 1 found without it.

## Test 4: Status vs homophily

- **Worlds:**

  | World | Education part of tastes |
  |---|---|
  | homophily | -1 x abs(edu difference) |
  | status | +2 x partner's education |
  | both | both terms |

- **Model:** well-specified + partner's education.
- **Claim:**
  - one-sided separates the two;
  - two-sided: the partner's-education coefficient shrinks toward zero,
    and the sorting that status competition produces is read as homophily.

## Outputs

- `results/layouts.pdf`
- `results/test1.csv` to `test4.csv`
- one figure per test
- `REPORT.md`: one table and one sentence per test.
