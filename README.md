# What a one-sided choice model sees in a two-sided partner market

Monte Carlo simulations accompanying the doctoral dissertation of Jesper
Lindmarker (Institute for Analytical Sociology, Linköping University). The
dissertation's papers estimate conditional logit models of partner choice
from Swedish register data. This repository asks what those models recover
when the partner market works the way partner markets do: both people have
to agree.

**Work in progress.** The October 2026 revision is described here. The
earlier simulations (March to September 2026) are in
`archive/2026_first_generation/`, with a README explaining what they did
and why they were replaced.

## The model, in one paragraph

A market of 1,000 men and 1,000 women in three ancestry groups, each person
with a continuous education score (and in some modules a location or an
institution). Everyone values a partner by a few pair characteristics (same
group, education difference, distance) plus a random taste drawn once and
held fixed. Partners are formed in one of three ways:

- **One-sided choice (`oracle`).** Everyone picks freely, and the chosen
  person has no say. This is exactly what conditional logit assumes, so it
  is the benchmark.
- **Gale-Shapley (`gs`).** Everyone meets everyone, both must prefer the
  match to staying single, and people compete for partners until no two
  would rather be together.
- **Sequential search (`seq`).** People who are still single meet a few new
  people each round, about 17 in all before pairing. The market clears among
  those who have met.

In both two-sided settings the value of staying single is calibrated so
that 85% pair up. The researcher then does what the papers do: takes each
partnered person, samples 99 alternatives from the opposite sex, and fits a
conditional logit.

## Modules

| Module | One question | Bears on |
|---|---|---|
| [A Recovery](modules/A_recovery) | Does the model recover its target under unequal groups, segregation and sampled choice sets? | All papers |
| [B Two-sided](modules/B_two_sided) | What can two-sided data not tell us: which sex holds a taste, and homophily vs status? | Papers 1, 4 |
| [C Mediation](modules/C_mediation) | Are KHB mediation shares trustworthy? | Papers 1, 4 |
| [D Decomposition](modules/D_decomposition) | Does the composition/coefficient split attribute change over cohorts correctly? | Paper 3 |

Each module folder has a `SPEC.md` (the claim, the design, the pass
criteria, written before the runs), the run and analysis scripts, a
`results/` folder, and a `REPORT.md`. The main results across modules are
summarised in [FINDINGS.md](FINDINGS.md).

## Layout

```
R/
  population.R   market: groups, education, space, institutions, utilities
  matching.R     oracle choice, deferred acceptance, sequential search, calibration
  estimate.R     choice sets (uniform, stratified + offset), clogit, KHB, predicted shares
  run_helpers.R  the shared settings (oracle, gs, gs_w, seq) and per-market code
tests/
  test_matching.R  stability, one-to-one, meeting budget, calibration (Rscript tests/test_matching.R)
modules/A_recovery ... modules/D_decomposition
archive/2026_first_generation/   superseded work, kept for the record
```

## Running

Requirements: R 4.x with `survival`, `parallel` and `splines` (base), and
`ggplot2` for figures. Run everything from the repository root:

```
Rscript tests/test_matching.R
Rscript modules/A_recovery/run_A.R
Rscript modules/A_recovery/analyse_A.R
```

and likewise for B to D.

Runs use a PSOCK cluster (separate R processes): forked workers crash
`survival::clogit` on macOS. Each cell has its own RNG stream, so results
do not depend on the number of workers. Every results file records the git
commit and seed it was produced with. On a 12-core laptop each module's
production run takes 15 to 40 minutes.

## Status

- A to D have been run and reported (October 2026). See `FINDINGS.md`.
- Not yet: a written-up paper, and a release with a DOI. Results may change
  as the work is reviewed.
