# Module A: Recovery

**Status:** spec revised 2026-10-13 after the mechanism validation. Production run: `run_A.R`.
Stylised toy model by design: no calibration to Swedish magnitudes.
**Replaces:** old Module 1 (group size), Module 2 (spatial), the oracle half of
Module 1b, and the Module 4 selectivity sweep.

## The one claim

> Conditional logit recovers its target under unequal group sizes, residential
> segregation and sampled or stratified choice sets. Under one-sided (oracle)
> choice the target is the chooser's own coefficient. Under two-sided matching
> the target is the pair surplus (the sum of both sides' coefficients), up to
> a scale factor that is the same for every coefficient.

A is also the control condition for B, C and D: if A passes, a distortion
found later comes from the departure that module introduces.

## The model

**Market.** N men and N women (default 1,000 each). Three ancestry groups
(default shares 70/20/10, same for both sexes). Education is one continuous
score, N(0, 1).

**Utility.** Man m's utility from woman w, and w's from m:

    U_mw = b_g * same_group + b_e * |edu_m - edu_w| + eps_mw
    V_wm = b_g * same_group + b_e * |edu_m - edu_w| + eta_wm

- Default coefficients: b_g = 1, b_e = -1.
- Both sexes have the same coefficients in A. Asymmetry between the sexes is
  Module B's subject.
- Shocks are i.i.d. Gumbel, drawn once per pair and held fixed.
- Optional outside option: u0 + Gumbel per person. The default is
  u0 = -Inf, meaning nobody prefers to stay single.

**Partnering: sequential search with deferred acceptance**
(`R/matching.R::search_market`; validated in `VALIDATION.md`).

- Tastes are fixed: utilities and outside options are drawn once.
- Each round, everyone still single meets new people of the other sex at
  random. Meetings are mutual and accumulate.
- After each round the market clears by deferred acceptance among everyone
  who has met. Couples stay provisional until the last round.
- The final matching is stable with respect to who met whom. Meeting
  everyone gives textbook Gale-Shapley.
- **The friction parameter is the absolute number of people met, not a share
  of N** (VALIDATION.md, V3).
- The outside option level u0 is calibrated per cell to **85% pairing**.
  That separates singlehood from friction.
- `tests/test_matching.R` (100+ checks) verifies:
  - one-to-one couples;
  - couples only among people who met;
  - no blocking pairs on the final network;
  - equivalence to Gale-Shapley at full budget;
  - the meeting budget;
  - the calibration.

The **oracle** (each ego picks freely, with no capacity limit and no say for
the other side) is kept as the reference, because it is the process clogit
assumes.

**Estimation (as in the papers).**

- One model per ego sex.
- Egos are partnered individuals only.
- Each choice set is the actual partner plus J - 1 = 99 alternatives drawn
  uniformly from the opposite-sex pool.
- Estimator: `survival::clogit`.
- With identical coefficients across sexes, the men's and women's estimates
  should agree. They do in the pilot.

## A0 pilot (2026-10-12, 20 reps per cell; one-shot meetings, superseded by VALIDATION.md)

`results/pilot_A0_summary.csv`. rho = estimate / target, where the target
is b for the oracle and 2b for deferred acceptance (DA).

| Mechanism | k | u0 | Pairing rate | rho, same group | rho, education |
|---|---|---|---|---|---|
| Oracle | n/a | none | 1.00 | 0.99 to 1.00 | 1.01 |
| DA, men propose | everyone | none | 1.00 | 0.94 | 0.93 |
| DA, women propose | everyone | none | 1.00 | 0.92 | 0.93 |
| DA, men propose, N = 500 / 2,000 | everyone | none | 1.00 | 0.87 / 0.94 | 0.88 / 0.95 |
| DA, men propose | 50 | none | 0.99 | 0.72 | 0.73 |
| DA, men propose | 10 | none | 0.95 | 0.42 | 0.44 |
| DA, men propose | 50 | 2 | 0.32 | 0.87 | 0.88 |
| DA, men propose | 10 | 2 | 0.11 | 0.86 to 0.90 | 0.87 to 0.88 |

What the pilot shows:

1. **The oracle recovers its target.** The pipeline works.
2. **Two-sided matching recovers the pair surplus up to a scale factor.**
   With everyone meeting everyone, rho is about 0.93 to 0.95 at
   N = 1,000 to 2,000. The proposing side barely matters.
3. **Search frictions shrink every coefficient by a similar factor.** Fewer
   meetings pull rho down (0.94, then 0.72, then 0.42), and within each cell
   rho is close for both terms (within 0.03 here; up to 0.07 in the fuller
   validation, so the shared scale is approximate).
   - Frictions act like extra noise in the logit scale: absolute levels fall
     and **ratios between coefficients survive**.
   - An outside option restores much of the level, because people can
     refuse a poor match instead of settling.
4. **The working target is therefore the pair surplus times an
   approximately common scale.** That is what A1 to A3 test the stability of.

Point 3 carries a warning for Module D. If the scale depends on market
frictions, a change in frictions across cohorts would look like a change in
every coefficient at once. D should test this directly.

## Blocks

| Block | Question | Varies |
|---|---|---|
| **A0 Anchor** | What does each mechanism recover with no threats? | Baseline: 70/20/10, no space, J = 100 uniform |
| **A1 Group size** | Does the common scale hold when one group dominates or is tiny? | Shares {33/33/33, 70/20/10, 85/10/5, 90/8/2} |
| **A2 Space** | Does ln(distance) absorb spatial confounding when groups are segregated? | Segregation {none, moderate, high} x model {with, without ln(dist)} |
| **A3 Sampling** | Do sampled and stratified choice sets recover the target? | J = 30 uniform; stratified (50 of 100 alternatives from ego's own group) with and without the McFadden-Manski offset. Separately: fastclogit vs `survival::clogit` on the stratified data |

**Settings, run in every block (A0 to A3):**

| Setting | Mechanism | Meetings | Outside option | Target |
|---|---|---|---|---|
| `oracle` | one-sided choice | everyone | off | b |
| `da_full` | deferred acceptance | everyone meets everyone (frictionless reference) | 85% pairing | 2b x scale |
| `da_seq` | sequential search, budget 50, 10 rounds (about 17 met; frictional default) | 50 | 85% pairing | 2b x scale |

**A2 adds space.** People sit on a unit square. Groups cluster around their
own centres, with the tightness set by the segregation level. Utility gains
b_d * ln(distance + 0.01), with b_d = -0.5, for both sides. Meetings stay
random with respect to distance: in A, space enters through tastes only.

**Replications:** 200 per cell.

## Pass criteria

| Mechanism | Passes if | Fails if |
|---|---|---|
| Oracle | abs(rho - 1) < 0.02 and 95% coverage in [0.93, 0.97] | Bias appears with small groups, segregation or stratified sampling |
| DA | Within each DA setting, each term's rho stays within 0.05 of its A0 value | rho shifts with group shares, segregation or sampling |

The gap in rho between terms is reported, not required to be zero: validation showed the shared scale is approximate (gap up to 0.07).

In A2, the model without ln(distance) is expected to fail under segregation.
That is a demonstration of omitted-variable bias, not a test of the
pipeline.

## Outputs

- `results/A{0..3}.csv`: one row per cell x coefficient x ego sex, carrying
  the commit and seed.
- `REPORT.md`: one pass/fail line per block, and one figure of rho across all
  cells.
- A comparison table against the old Modules 1, 2 and 4.

## Not in A

- Sexes with different coefficients and the direction of assortativity (B).
- Homophily versus status (B).
- Mediation (C).
- Change over time (D).
