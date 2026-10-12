# Module A: Recovery

**Status:** spec revised 2026-10-12 after Jesper's comments and the A0 pilot.
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

**Partnering: deferred acceptance among the people you meet**
(`R/matching.R`).

- Each person meets k random people of the other sex. Meetings are mutual and
  unrelated to traits.
- Proposers work down their list of the people they met, best first, and
  skip anyone below their outside option.
- Receivers hold the best proposal so far and drop it if a better one comes.
- The process runs until no proposer has anyone left to try.
- The result is stable with respect to the meeting network. Setting
  k = everyone gives textbook Gale-Shapley.
- `tests/test_matching.R` verifies: one-to-one matching, couples only among
  people who met, no blocking pairs, and an identical result from identical
  utilities.

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

## A0 pilot (2026-10-12, 20 reps per cell)

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
3. **Search frictions shrink every coefficient by the same factor.** Fewer
   meetings pull rho down (0.94, then 0.72, then 0.42), but within each cell
   rho is the same for both terms, to within 0.03.
   - Frictions act like extra noise in the logit scale: absolute levels fall
     and **ratios between coefficients survive**.
   - An outside option restores much of the level, because people can
     refuse a poor match instead of settling.
4. **The working target is therefore the pair surplus times a common
   scale.** That is what A1 to A3 test the stability of.

Point 3 carries a warning for Module D. If the scale depends on market
frictions, a change in frictions across cohorts would look like a change in
every coefficient at once. D should test this directly.

## Blocks

| Block | Question | Varies |
|---|---|---|
| **A0 Anchor** | What does each mechanism recover with no threats? | Done as pilot; production run at 200 reps |
| **A1 Group size** | Does the common scale hold when one group dominates or is tiny? | Shares {33/33/33, 70/20/10, 85/10/5, 90/8/2} |
| **A2 Space** | Does ln(distance) absorb spatial confounding when groups are segregated? | Segregation {none, moderate, high} x model {with, without ln(dist)} |
| **A3 Sampling** | Do sampled and stratified choice sets recover the target? | J {30, 100}; stratified (oversampled same-group alternatives) with and without the McFadden-Manski offset; fastclogit vs `survival::clogit` on the stratified data |

**Mechanisms in A1 to A3:**

- Oracle.
- DA at two friction levels:
  - k = everyone (frictionless);
  - k = 50 (frictional).
- Outside option off.

Running both DA levels tests whether the threats interact with frictions.

**A2 adds space.** People sit on a unit square. Groups cluster around their
own centres, with the tightness set by the segregation level. Utility gains
b_d * ln(distance), with b_d = -0.5, and both sides get the same term, as
above.

**Replications:** 200 per cell.

## Pass criteria

| Mechanism | Passes if | Fails if |
|---|---|---|
| Oracle | abs(rho - 1) < 0.02 and 95% coverage in [0.93, 0.97] | Bias appears with small groups, segregation or stratified sampling |
| DA | Within each DA setting, rho stays within 0.05 of its A0 value **and** is equal across terms (within 0.05) | rho shifts with group shares, segregation or sampling, or differs between terms |

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
