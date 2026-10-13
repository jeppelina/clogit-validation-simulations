# Findings

**October 2026 revision.** Detail, numbers and robustness checks are in each
module's `REPORT.md`.

## The setup

- A stylised market: 1,000 men and 1,000 women, two groups, continuous
  education, fixed tastes.
- Partners are formed either one-sidedly, which is what conditional logit
  assumes, or two-sidedly (Gale-Shapley: both must agree).
- The model is then fitted the way the dissertation's papers fit it.
- The one-sided runs recover the true tastes in every module. That checks
  the code. Everything below is about what two-sidedness does.

## Five claims

**1. A well-specified conditional logit does what it aims to do, up to a
scale set by search frictions** (Module A, test 1 and 2).
- Well-specified means:
  - egos from every group pooled;
  - alter-group constants;
  - distance in the model when there is space.
- Under one-sided choice it recovers the true tastes.
- Under two-sided matching it recovers the pair-level assortativity (what
  a log-linear model of the couple table measures) at a stable scale:
  about 0.96 when everyone meets everyone, lower with search frictions.
- It is free of group size (50% to 5% minority) and of residential
  segregation (two layouts, moderate and strong). Without distance,
  propinquity is read as endogamy: up to +25%.
- **Absolute coefficient sizes say as much about the market as about
  tastes.**

**2. It cannot say which sex holds a taste** (Module A, test 3).
- Worlds where only men, only women, or both care give identical men's and
  women's estimates.
- With the alter-group constant in the model, an unbalanced sex ratio does
  not create a spurious sex gap. It does without the constant.

**3. Status preferences are nearly invisible, and what they do produce
looks like homophily** (Module A, test 4).
- A status taste of 2 is estimated at 0.22.
- The educational sorting that status competition produces (correlation
  0.35) is reported as homophily: -0.56 where the truth is 0.
- **Educational homogamy estimated from couples is evidence of sorting,
  not of a taste for similarity.**

**4. Over time, the well-specified model separates assortativity from
composition** (Module D).
- The minority grows from 10% to 30% across two cohorts.
- The same-group part of an Oaxaca-type decomposition is about 0 when
  tastes are constant, and within 7% of the true effect when they harden or
  soften, under one-sided and two-sided matching.
- The alter-group constant absorbs the two-sided market's own adjustment to
  composition.
- Fitted from one group's egos only, as in Paper 3, the same decomposition
  reports +3 to +4 points of spurious hardening (Module D v1). **Use the
  whole couple table.**

**5. Unobserved confounding in a mediation analysis lands mainly on the
setting coefficient** (Module C).
- An unobserved trait sorts people into settings and also shapes partner
  choice. The setting's coefficient is inflated by about 0.3 to 0.55 (in
  units of the homogamy coefficient). A setting with no effect gets 0.56.
- The direct homogamy coefficient moves by at most 10%.
- The mediation share is inflated by a few points. That is large when the
  real channel is weak (0.13 vs 0.08), modest when it is strong (0.23 vs
  0.19).
- Two-sided matching does not change where the bias lands.

## What this does and does not show

- It is a toy. It shows that these mechanisms exist and roughly how large
  they can be. It does not show how large they are in Sweden.
- The two-sided rule is Gale-Shapley with everyone meeting everyone, plus
  an outside option (85% pair). Search frictions lower the scale of every
  coefficient (Module A validation); they are not varied in the main tests.
- Space and settings enter through tastes, not through who meets whom.
- The papers' main descriptive conclusions survive:
  - composition drives the rise in exogamy;
  - propinquity explains a minority of endogamy;
  - education sorting is strong.

  What needs care:
  - reading coefficients as individual tastes;
  - sex differences;
  - status vs homophily;
  - setting coefficients under selection into settings;
  - any comparison across groups or over time made without alter-group
    constants on the whole couple table.
