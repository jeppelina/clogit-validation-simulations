> **Superseded (2026-10-14).** Module A version 3 now carries these tests
> (test 3: sex differences; test 4: status vs homophily) with the
> well-specified model. Module A v3 shows that this report's sex-ratio gap
> (B1b) disappears once an alter-group constant is in the model. Kept for
> the record.

# Module B report: what two-sided data cannot tell us

**Run:** 2026-10-13. `run_B.R` (200 reps per cell), `analyse_B.R`,
`diag_B2.R` (100 reps per cell). Settings:
- `oracle`: one-sided choice;
- `gs`: Gale-Shapley, men proposing;
- `gs_w`: Gale-Shapley, women proposing;
- `seq`: sequential search, budget 50.

Every two-sided setting uses 85% pairing.

## B1a. Which sex holds the assortativity? The data cannot say

The pair surplus on same group is 2 in every world. Estimated same-group
coefficient:

| World | Setting | Men's model | Women's model |
|---|---|---|---|
| Only men care (2, 0) | oracle | **1.99** | **0.00** |
| Only women care (0, 2) | oracle | **0.01** | **2.01** |
| Both care (1, 1) | oracle | 0.99 | 1.01 |
| Only men care | gs | 1.81 | 1.80 |
| Only women care | gs | 1.81 | 1.82 |
| Both care | gs | 1.90 | 1.90 |
| Only men care | gs_w | 1.79 | 1.79 |
| Only women care | gs_w | 1.81 | 1.81 |
| Only men care | seq | 0.91 | 0.91 |
| Only women care | seq | 0.91 | 0.90 |
| Both care | seq | 1.23 | 1.23 |

Monte Carlo SEs are 0.005 to 0.008.

- **One-sided choice recovers each side's own taste.**
- **Under every two-sided setting, the men's and women's estimates are
  identical** (differences of 0.006 or less, within simulation noise), and
  identical between "only men" and "only women".
  - Which sex proposes makes no difference.
  - A one-sided clogit fitted to couples estimates a property of the pair,
    not of the side whose choices it models.
- **The level depends on how the surplus is split, not only on its size.**
  A taste held by one side yields a smaller coefficient than the same
  surplus shared by both: 1.81 against 1.90 under `gs`, and 0.91 against
  1.23 under `seq`. With search frictions, the side that does not care
  accepts whoever it meets first.

## B1b. Identical tastes, different pools: sex differences appear anyway

Both sexes have identical tastes (same group 1, education -1). Same-group
coefficient:

| World | Setting | Men | Women | Difference |
|---|---|---|---|---|
| symmetric pools | gs | 1.90 | 1.90 | -0.00 |
| minority groups 60% male | oracle | 1.00 | 1.00 | -0.00 |
| minority groups 60% male | gs | 1.97 | 1.58 | **+0.39** (MC SE 0.011) |
| minority groups 60% male | seq | 1.33 | 1.05 | **+0.28** (MC SE 0.009) |
| women 0.5 SD more educated | gs | 1.89 | 1.89 | +0.00 |

- **An unbalanced sex ratio within groups produces a large sex difference
  in the same-group coefficient**, with identical tastes. In this world
  the men's coefficient is 25% larger than the women's.
- The oracle shows no difference, so this comes from two-sided matching.
  It is the same group-size mechanism as Module A: within the minority
  groups, men and women face different relative pool sizes.
- **The education gap did not produce a sex difference** on the continuous
  education-difference term, under any setting.
  - Module 8's earlier finding of a sex gap under competition involved a
    categorical top-education term with rationing at the top. It ran on the
    old, flawed matching algorithm.
  - This run does not support a general claim that education composition
    creates sex gaps.

## B2. Status preferences are hidden, and what survives looks like homophily

Estimated model: same group + abs(edu difference) + partner's education.
Men's models (women's are symmetric):

| World (true tastes) | Setting | abs(edu diff) | Partner's edu | Partners' edu correlation |
|---|---|---|---|---|
| Homophily (-1, 0) | oracle | -1.00 | 0.00 | 0.52 |
| Homophily | gs | -1.93 | 0.00 | 0.78 |
| Status (0, 1) | oracle | 0.00 | **1.00** | 0.00 |
| Status | gs | -0.14 | **0.14** | 0.09 |
| Status | seq | -0.20 | 0.13 | 0.13 |
| Mixed (-1, 1) | oracle | -1.00 | 1.00 | 0.44 |
| Mixed | gs | -1.93 | **0.17** | 0.77 |
| Strong status (0, 3) | oracle* | -0.55 | 3.55 | 0.00 |
| Strong status | gs | **-1.16** | 0.24 | **0.61** |

\* With a status taste of 3, many men choose the same few women under
oracle. The fit is close to separation and is itself off (3.55 against 3).
Treat magnitudes at that strength as rough.

1. **Status preferences largely disappear.** The coefficient on partner's
   education drops from 1.00 to 0.14 (status world) and to 0.17 (mixed
   world), against a true 1.
   - In the mixed world, homophily is recovered as fully as in the pure
     homophily world (-1.93), while the status taste is nearly invisible.
     A researcher would conclude "homophily, little status".
2. **Why.** When everyone wants the same partners, matching is decided by
   competition, and competition is not visible in a comparison of partners
   with random alternatives.
   - A status taste adds the same amount to every pair containing a given
     person, whoever the partner is. So it mostly decides who gets matched
     at all and who stays single. Partnered men are 0.9 SD more educated
     than single men in the status world.
   - It does not decide who is matched with whom.
3. **When the status taste is strong enough to sort couples,** competition
   produces real educational sorting (correlation 0.61), and the model
   reports it as homophily: -1.16 against a true 0.

**Claim, stated carefully.** A one-sided model fitted to two-sided data
cannot recover status preferences. Whatever sorting they produce through
competition is read as homophily. The converse does not hold: true homophily
is not mistaken for status (partner's education is 0.00 in the homophily
world).

## For the papers

- **Per-sex estimates (Papers 1 and 4).** A difference between the sexes'
  coefficients cannot be read as a difference in tastes. It can arise from
  sex-ratio differences within groups alone (B1b), and identical estimates
  are consistent with only one sex caring (B1a).
- **Homogamy as homophily (Paper 4).** Educational homogamy estimated from
  couples is consistent with status competition. It is evidence of sorting,
  not of a taste for similarity.
