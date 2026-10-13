# What does a conditional logit of partner choice estimate when both partners have to agree?

Monte Carlo simulations of conditional logit models of partner choice, the
workhorse of register-based research on assortative mating. The model treats
partner choice as one-sided: each person picks a partner from a set of
alternatives. Real unions are two-sided: both people have to agree. These
simulations generate partner markets with known tastes, form unions in both
ways, and check what the model recovers.

## TL;DR

**A. A well-specified model does what it should, up to a point.**
- A conditional logit with **alter-group constants** recovers same-group
  assortativity regardless of **group size**, and regardless of
  **residential segregation** as long as **ln(distance)** is in the model.
- Under two-sided matching it estimates the **pair's** assortativity, both
  partners' tastes combined, at a stable scale. That is the quantity a
  log-linear model of the couple table measures.
- **It cannot:**
  - say **which sex** holds a taste: models for men and for women give the
    same estimate whether only men, only women, or both care;
  - see **status preferences**: a taste for a more educated partner shows
    up at about a tenth of its size, and the sorting it produces is read as
    homophily.
- **Without the alter-group constant**, small groups' assortativity is
  underestimated by about 15% for a 5% minority. **Without distance**, propinquity
  is read as endogamy (up to +26%).

**B. In a mediation analysis, unobserved confounding lands on the setting.**
- If an unobserved trait drives both who ends up in a setting (a workplace,
  say) and whom people choose, the **setting's coefficient** is inflated. A
  setting with no effect gets a clearly positive one.
- The **direct homogamy** coefficient moves by at most about 10%.
- The **mediation share** is inflated by 2 to 5 points: a large error when
  the setting's real role is small.

**C. Over time, the well-specified model separates assortativity from
composition.**
- As a minority grows from 10% to 30%, the part of a decomposition
  attributed to assortativity is about zero when tastes are constant.
- It tracks real hardening or softening within 8%.
- The alter-group constant absorbs the market's own adjustment to the new
  composition.

## The setup in brief

- **Market:**
  - 1,000 men and 1,000 women;
  - a majority and a minority group;
  - education ~ N(0, 1);
  - in some tests, locations, settings, or an unobserved trait.
- **Tastes:** a linear utility in pair characteristics (same group,
  |education difference|, ...) plus a random taste for each pair, drawn
  once.
- **Unions, formed in two ways:**
  - **One-sided choice:** each person picks their best option. This is the
    process the conditional logit assumes, and the benchmark.
  - **Two-sided matching (Gale-Shapley):** men propose, women hold the best
    offer, and nobody accepts a partner worse than staying single. The
    value of staying single is set so that 85% pair up.
- **Estimation, as in register studies:**
  - each partnered person's choice set is their partner plus 99 people of
    the other sex drawn at random;
  - one conditional logit with egos from all groups;
  - the **well-specified** model has a same-group term, a constant for the
    alternative being from the minority, and |education difference| (plus
    terms specific to each study).
- 200 simulated markets per cell.

## Repository

```
R/                    core library
  population.R        markets, layouts, settings, utilities
  matching.R          one-sided choice, Gale-Shapley, calibration of pairing
  estimate.R          choice sets, conditional logit, KHB, predicted shares
  run_helpers.R       the two settings, forming unions, the parallel cluster
tests/test_core.R     tests of the core library
A_recovery/           group size, segregation, sex differences, status vs homophily
B_mediation/          mediation under unobserved confounding
C_change_over_time/   decomposing change across two cohorts
```

Each study folder has a README with its design, results and figures, a
`run.R` that simulates, an `analyse.R` that summarises, and the outputs in
`results/`.

## Running

R 4.x with `survival` and `ggplot2` (and the base packages `parallel` and
`stats`). From the repository root:

```
Rscript tests/test_core.R
Rscript A_recovery/run.R && Rscript A_recovery/analyse.R
Rscript B_mediation/run.R && Rscript B_mediation/analyse.R
Rscript C_change_over_time/run.R && Rscript C_change_over_time/analyse.R
```

- **Parallelism.** Simulations run in parallel on separate R processes,
  using all but two cores. Each cell has its own random number stream,
  given by the `seed` column in the results.
- **Run time.** About 15, 12 and 9 minutes for A, B and C on a 12-core
  laptop.

## Citation

Lindmarker, Jesper. 2026. *What does a conditional logit of partner choice
estimate when both partners have to agree? Simulation code and results.*
Institute for Analytical Sociology, Linköping University.
