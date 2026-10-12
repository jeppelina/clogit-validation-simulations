# Module B: What two-sided data cannot tell us

**Status:** spec, 2026-10-13. Uses the market, partnering rule and settings
validated in Module A.

## The two claims

> **B1. Direction.** In a two-sided market, a one-sided conditional logit
> cannot say which sex holds the assortativity. Worlds in which only men,
> only women, or both care about ancestry give the same per-sex estimates,
> as long as the pair surplus is the same. Conversely, men and women with
> identical tastes can get different per-sex estimates when the two sexes'
> pools differ in composition.
>
> **B2. Homophily versus status.** When everyone prefers a more educated
> partner and nobody cares about similarity, two-sided competition still
> sorts people by education. The model then reports homophily.

Both are about the meaning of a coefficient, not its accuracy. Module A
showed the model recovers the pair surplus (up to scale); Module B shows
what the pair surplus hides.

## Shared setup

Same as Module A:
- 1,000 men and 1,000 women;
- three groups, 70/20/10;
- continuous education;
- one model per ego sex, 100 alternatives drawn uniformly.

Settings:

| Setting | What it is |
|---|---|
| `oracle` | One-sided choice. The reference: it recovers each side's own tastes. |
| `gs` | Textbook Gale-Shapley, men proposing (everyone meets everyone), outside option calibrated to 85% pairing. The main two-sided setting. |
| `gs_w` | As `gs`, women proposing. Only where direction is the question (B1). |
| `seq` | Sequential search, budget 50, 10 rounds, 85% pairing. A robustness check on the main cells. |

200 replications per cell.

## B1. Direction

**B1a. Who holds the assortativity?** Same-group coefficient by sex. The
pair surplus is 2 in every world:

| World | Men's b_g | Women's b_g |
|---|---|---|
| men_only | 2 | 0 |
| women_only | 0 | 2 |
| both | 1 | 1 |

Education coefficient -1 for both sexes in every world.

- **Expected under oracle:** per-sex estimates recover each side's own
  coefficient (2/0, 0/2, 1/1).
- **Claim under `gs`:** both sexes' estimates are about equal and the same
  across the three worlds.
- **What would falsify it:** per-sex estimates under `gs` that track their
  own side's coefficient.

**B1b. Same tastes, different pools.** Both sexes have identical tastes
(b_g = 1, b_e = -1). Two compositional asymmetries, each against a symmetric
baseline:

- **sex_ratio:** the minority groups are 60% male and the majority balances
  the totals, as in many migrant-origin groups.
- **edu_gap:** women's education is shifted up by 0.5 SD, as in Sweden.

- **Claim:** under `gs`, men's and women's estimates diverge on the term
  where the pools differ: same_group for sex_ratio, edu_diff for edu_gap.
  Under oracle they don't.
- **What would falsify it:** no divergence beyond simulation noise. Then the
  claim narrows to B1a alone.

## B2. Homophily versus status

Three worlds, same for both sexes, plus b_g = 1:

| World | Utility (education part) |
|---|---|
| homophily | -1 x abs(edu_m - edu_w) |
| status | +1 x partner's education |
| mixed | both terms |

Estimated model: same_group + abs(edu difference) + partner's education.

- **Expected under oracle:** each world is recovered. The absolute
  difference and the partner's level are separately identified, as in the
  old Module 1b.
- **Claim under `gs`:** in the status world, competition sorts people by
  education, so the estimate on abs(edu difference) is clearly negative,
  i.e. it looks like homophily. The coefficient on partner's education is
  attenuated (with 85% pairing, the least educated stay single, so some
  signal survives).
- **What would falsify it:** the status world under `gs` shows no
  homophily-like coefficient.

## Outputs

- `results/B1a.csv`, `B1b.csv`, `B2.csv`;
- `REPORT.md` with one table per claim;
- one figure: per-sex estimates across the B1a worlds.
