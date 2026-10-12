# Validation of the partnering mechanism (2026-10-13)

Script: `validate_search.R`. Results: `results/validate_search_summary.csv`
(raw: `validate_search_raw.csv`). N = 1,000 per sex unless stated; 20 reps
per cell (50 for the null). rho = estimate / pair surplus (2b) for
deferred acceptance, estimate / b for the oracle. "Met" is the mean
number of people a man actually met before the search ended.

## V1 Null: no sorting is manufactured

With all coefficients at zero, the mean estimates are 0.001 to 0.012 under
every mechanism (oracle; everyone meets everyone; sequential search with
85% pairing). Coverage of 0 is 0.92 to 1.00. Majority endogamy is 0.70,
which is the majority's share of the market. **Pass.**

## V2 Grid: rho is set by how many people a person actually meets

| Search | Budget | Outside option | Met | Pairing | rho, group | rho, edu |
|---|---|---|---|---|---|---|
| Sequential (10 rounds) | 10 | off | 2.5 | 0.94 | 0.17 | 0.18 |
| Sequential | 10 | 85% | 3.7 | 0.85 | 0.30 | 0.35 |
| Sequential | 50 | off | 6.3 | 0.99 | 0.32 | 0.35 |
| One-shot | 10 | off / 85% | 10 | 0.95 / 0.86 | 0.43 / 0.49 | 0.44 / 0.54 |
| Sequential | 50 | 85% | 17 | 0.84 | 0.61 | 0.67 |
| Sequential | 200 | off | 21 | 1.00 | 0.55 | 0.59 |
| One-shot | 50 | off / 85% | 50 | 0.99 / 0.86 | 0.71 / 0.75 | 0.72 / 0.79 |
| Sequential | 200 | 85% | 64 | 0.85 | 0.81 | 0.83 |
| One-shot | 200 | off / 85% | 200 | 1.00 / 0.85 | 0.86 / 0.88 | 0.86 / 0.90 |
| Everyone meets everyone | 1,000 | off / 85% | 1,000 | 1.00 / 0.85 | 0.90 / 0.93 | 0.91 / 0.98 |
| Oracle | n/a | off | n/a | 1.00 | 1.02 | 0.99 |

1. **rho rises with the number of people actually met:** about 0.2 at 3
   meetings, 0.6 at 17, 0.75 at 50, 0.9 at everyone.
   - Budget, search format and outside option matter almost only through
     that number.
   - At a given number met, the outside option adds a little (0.03 to 0.06).
2. **Sequential search meets far fewer people than its budget.** People stop
   searching once paired.
   - Without an outside option, anyone acceptable is accepted, so people
     pair with one of their first few contacts (6 meetings on a budget of 50).
   - This is the "settle early" market, and it gives the strongest
     attenuation.
3. **The two coefficients shrink by similar but not identical factors.**
   - The education term keeps a little more of its size than the group term:
     a gap of 0.01 to 0.06 in most cells, up to 0.07 at N = 500.
   - So the "common scale" from the pilot is approximate, not exact.
   - A plausible reason: with few meetings, minority members rarely meet
     anyone from their own group, which costs the group term more than the
     continuous education term.

## V3 Scaling: the absolute number met matters, not the share of N

Sequential search, 85% pairing:

| N | Budget | Share of N | Met | rho, group |
|---|---|---|---|---|
| 500 | 50 | 10% | 16.9 | 0.58 to 0.59 |
| 1,000 | 50 | 5% | 17.0 | 0.61 to 0.62 |
| 2,000 | 50 | 2.5% | 16.5 | 0.59 to 0.60 |
| 500 | 25 | 5% | 8.7 | 0.49 |
| 2,000 | 100 | 5% | 32.9 | 0.70 to 0.71 |

- A fixed budget gives the same rho at every N.
- A fixed share does not: rho rises with N because more people are met.
- **So the friction parameter is the absolute number of people met**, which
  is also the more realistic reading: a person meets a limited number of
  potential partners, whatever the size of the country.

## What this settles for Modules A to D

- The mechanism manufactures no sorting (V1).
- Two-sided matching estimates the pair surplus times an attenuation factor
  set mainly by search frictions. The factor is shared closely, but not
  exactly, across coefficients.
- Coverage against the unscaled pair surplus is near 0 in every two-sided
  cell. That is expected given rho < 1, and is not a pipeline failure. The
  meaningful test in A1 to A3 is whether rho stays where A0 puts it.
