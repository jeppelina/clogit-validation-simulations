# Module A report: recovery

**Run:** 2026-10-13. `run_A.R` (200 reps per cell), `analyse_A.R`,
`diag_A.R` (40 reps per cell), `check_fastclogit.R`. Settings as in
`SPEC.md`:
- `oracle`: one-sided choice;
- `da_full`: everyone meets everyone, 85% pairing;
- `da_seq`: sequential search, budget 50, 85% pairing.

rho = estimate / target. The target is b under oracle and the pair surplus
2b under two-sided matching. Both ego sexes are pooled; their estimates
agree to within 0.02 in every cell.

## Verdict

| Claim | Oracle | Two-sided |
|---|---|---|
| Recovers its target with balanced groups | Pass (rho 1.00) | Pass, up to a scale: 0.96 (`da_full`), 0.61 to 0.66 (`da_seq`) |
| Sampling designs (J = 30; stratified with the McFadden-Manski offset) | Pass | Pass (within 0.02 of A0) |
| Unequal group sizes | Pass | **Fail for the same-group term** |
| Segregation, distance controlled | Pass | **Fail for the same-group and distance terms** |
| fastclogit equals `survival::clogit` | Max coefficient difference 2e-7 over 120 fits | |

The oracle passes everywhere. So does the pipeline: estimation, sampling,
the offset and the kernel. The failures are properties of two-sided
matching, not of the code, and they are the module's main result.

## 1. With balanced groups, two-sided matching estimates the pair surplus times a scale

| Setting | rho, same group | rho, education |
|---|---|---|
| oracle | 1.003 | 1.003 |
| da_full | 0.961 | 0.961 |
| da_seq | 0.607 | 0.660 |

The scale is set by search frictions: how many people a person meets
(`VALIDATION.md`). It is shared closely, but not exactly, across terms.

**Sampling.** Results relative to A0:
- J = 30: within 0.02.
- Stratified sampling with the offset: within 0.01.
- Stratified sampling **without** the offset: the same-group term is biased
  down by 0.12 under oracle (0.884) and by 0.14 under `da_full`. The offset
  is necessary and sufficient.

## 2. Under two-sided matching, the same-group coefficient depends on group sizes

Pooled same-group rho as the majority grows:

| Shares | oracle | da_full | da_seq |
|---|---|---|---|
| 33/33/33 | 1.003 | 0.960 | 0.641 |
| 70/20/10 | 1.003 | 0.961 | 0.607 |
| 85/10/5 | 1.003 | 0.901 | 0.536 |
| 90/8/2 | 1.000 | 0.860 | 0.506 |

The education term does not move (0.957 to 0.969 under `da_full`).

`diag_A.R` splits the coefficient by the ego's group. Tastes are identical
for everyone, yet:

| Shares | Setting | Majority (g1) | Middle group (g2) | Smallest group (g3) |
|---|---|---|---|---|
| 70/20/10 | gs | 0.70 | 1.13 | 1.26 |
| 90/8/2, N = 1,000 | gs | 0.49 | 1.59 | not stable (separation in some reps) |
| 90/8/2, N = 3,000 | gs | 0.49 | 1.55 | not stable |
| 90/8/2, N = 1,000 | seq | 0.34 | 0.99 | not stable |

**Larger groups look less endogamous and smaller groups look more
endogamous**, by a wide margin, with no difference in tastes.

- The effect is the same at N = 3,000 as at 1,000, so it is not a
  small-numbers artefact.
- The pooled coefficient falls as the majority grows because the majority
  dominates the pool of egos.

Why: conditional logit takes the population as the opportunity set. It
expects a group's endogamy to scale with the group's share of the pool. In
a two-sided market, members of a small group who want each other can
largely find each other, because the wanting is mutual. A large group's
members are spread over many partners who also have other options.
Observed endogamy therefore responds to group size less than the logit
assumes, and the coefficient absorbs the difference.

**Consequences:**
- **Comparisons of the same-group coefficient across groups of different
  size** (Papers 1 and 2) confound tastes with group size.
- **When composition changes over time** (Paper 3), the coefficient changes
  even if tastes do not. Module D tests this directly.

### Why, in one table (`explain_group_size.R`, 20 reps, two groups, everyone pairs)

| Shares | Setting | Majority's coefficient | Minority's coefficient | Sum | Couple-table log odds ratio | Majority endogamy | Minority endogamy |
|---|---|---|---|---|---|---|---|
| 50/50 | oracle | 0.98 | 0.98 | 1.96 | 1.97 | 0.73 | 0.73 |
| 70/30 | oracle | 1.01 | 0.99 | 2.00 | 2.00 | 0.87 | 0.54 |
| 90/10 | oracle | 1.00 | 1.00 | 2.01 | 2.00 | 0.96 | 0.23 |
| 50/50 | gs | 1.89 | 1.90 | 3.79 | 3.79 | 0.87 | 0.87 |
| 70/30 | gs | 1.53 | 2.24 | 3.76 | 3.76 | 0.91 | 0.80 |
| 90/10 | gs | 0.95 | 2.72 | 3.67 | 3.66 | 0.96 | 0.63 |

**Under oracle**, each group's coefficient is about 1 at every composition.
The logit's adjustment for group size is exactly right when choices are
one-sided.

**Under Gale-Shapley:**
- Mixed couples are shared: each one is an out-marriage for one majority
  person and one minority person. In the table the two kinds of mixed
  couple are always equal in number (37.5 and 37.5 at 90/10). Under oracle
  they are not (36 and 77).
- The same 37 mixed couples are 4% of majority men's unions but 37% of
  minority men's.
- The logit judges each group's endogamy against that group's own share of
  the pool. So the minority (63% endogamous against 10% availability) gets
  a large coefficient, and the majority (96% against 90%) a small one.
  Tastes are the same.
- **The sum of the two coefficients equals the couple-table log odds ratio
  and is stable across compositions** (3.79, 3.76, 3.67). Group size does
  not change how strongly the pair is assortative. It changes how that is
  split between the groups' coefficients.

## 3. Under two-sided matching and segregation, distance does not fully absorb spatial endogamy

Distance is in the model (`with_dist`):

| Segregation | Setting | rho, same group | rho, ln(distance) |
|---|---|---|---|
| 0 | da_full | 0.939 | 0.922 |
| 0.9 | da_full | 1.032 | 0.820 |
| 0 | da_seq | 0.590 | 0.566 |
| 0.9 | da_seq | 0.735 | 0.475 |
| 0.9 | oracle | 0.998 | 1.009 |

- As segregation rises, the distance coefficient loses scale and the
  same-group coefficient gains it.
- Replacing linear ln(distance) with a 5-df spline does not change this
  (`diag_A.R`: 1.036 to 1.068 under gs; 0.743 in both specifications
  under seq), so it is not a functional-form problem.
- A plausible reading, which these runs have not tested: segregation gives
  each group a local market where it is the majority, and Section 2's
  group-size effect operates on that local composition.

**Consequence:** under two-sided matching, residual endogamy net of distance
is overstated, and the share attributed to propinquity is understated
(quantified in Module C, cell C4).

## Not covered

- Paper 4's widest offsets: the kernel check reached a within-set spread of
  5.1, against 10.5 in Paper 4.
- The proposing side in A1 to A3. The pilot showed it barely matters with
  balanced groups.
