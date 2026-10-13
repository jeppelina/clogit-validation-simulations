# Findings

**October 2026 revision.** Detail, numbers and robustness checks are in each
module's `REPORT.md`.

## The setup

- A stylised market: 1,000 men and 1,000 women, three groups, continuous
  education, fixed tastes.
- Partners are formed either one-sidedly, which is what conditional logit
  assumes, or two-sidedly (Gale-Shapley, or sequential search with
  deferred acceptance).
- The model is then fitted the way the dissertation's papers fit it.
- The one-sided runs recover the true tastes in every module. That checks
  the code. Everything below is about what two-sidedness does.

## Five claims

**1. A well-specified conditional logit does what it aims to do, up to a
scale set by search frictions** (Module A).
- Well-specified means:
  - egos from every group pooled;
  - one same-group term per group;
  - alter-group constants.
- Under one-sided choice it recovers the true tastes.
- Under two-sided matching it recovers the same quantity as a log-linear
  model of the couple table: the pair-level assortativity, equal across
  groups and stable from 33/33/33 to 85/10/5. Choice-set size and
  stratified sampling (with the offset) do not matter.
- The scale depends on how many potential partners people actually meet:
  0.95 when everyone meets everyone, about 0.65 at about 17 meetings.
  **Absolute coefficient sizes say as much about the market as about
  tastes.**
- **Limit:** under heavy residential segregation, two-sided matching shifts
  about 10% of the weight from distance to the same-group terms, even in
  this specification.

**2. It cannot say which sex holds a taste** (Module B).
- Worlds where only men, only women, or both care about ancestry give
  identical men's and women's estimates.
- With identical tastes, an unbalanced sex ratio within groups produces a
  25% sex difference in the same-group coefficient.
- **Sex differences in coefficients are not evidence of sex differences in
  tastes.**

**3. Status preferences are invisible, and what they do produce looks like
homophily** (Module B).
- A taste for a more educated partner shrinks to about 15% of its value.
  In a world with both homophily and status, only homophily is seen.
- When status tastes are strong enough to sort couples, the model reports
  homophily.
- **Educational homogamy estimated from couples is evidence of sorting,
  not of a taste for similarity.**

**4. A same-group coefficient estimated from one side of the market depends
on group size, so a shrinking majority looks like hardening** (Modules A and
D).
- With identical tastes and no alter-group constants, larger groups get
  smaller same-group coefficients and smaller groups larger ones. At 90/8/2:
  majority 0.49, 8% group 1.59, relative to the pair surplus.
- **This is a specification issue, not a flaw of conditional logit.**
  - Pooling egos from all groups and adding alter-group constants (the
    column margin) removes it.
  - That model gives the same numbers as a log-linear model of the couple
    table, whose odds ratios are free of group size, as the log-linear
    tradition says.
- A design that uses egos from one group only cannot include that margin.
- In such a design (Paper 3's: Swedish-origin egos), when the majority
  shrinks across cohorts the composition/coefficient decomposition reports
  a coefficient change toward endogamy (+3 to +4 points of majority
  endogamy) even when tastes do not change. That is 23 to 51% of the
  composition component, depending on taste strength.
- The effect is robust to:
  - which sex proposes;
  - the outside option;
  - market size;
  - which sex is the ego;
  - Paper 3's alter-category specification.
- Real changes in tastes keep their sign when they are larger than the
  bias. The bias is an upward shift, so a real softening smaller than it
  would be reported as hardening.
- An independent code review re-implemented the pipeline and reproduced
  the result. It persists with no sampling of alternatives.
- **The couple-table log odds ratio (or a log-linear association term)
  stays flat when tastes are flat and moves when they move.** It is a check
  that can be run on the register data.

**5. KHB mediation shares are clean with a well-specified homogamy term;
coarse measures and two-sided matching distort them** (Module C).
- A channel that does nothing gets a share of zero when education enters
  the model as it enters tastes.
- With a coarse homogamy measure ("same tercile"), the same channel gets
  8%, and a real channel's share doubles (0.22 to 0.52, and 0.71 under
  Gale-Shapley).
- Two-sided matching understates propinquity's share of endogamy: 0.37
  becomes 0.31 to 0.27 under heavy segregation.

## What this does and does not show

- It is a toy. It shows that these mechanisms exist and roughly how large
  they can be. It does not show how large they are in Sweden.
- The partnering rule is stable matching among the people one meets.
  Meetings are random with respect to space and institutions, so space and
  institutions enter through tastes only.
- The papers' main descriptive conclusions survive:
  - composition drives the rise in exogamy;
  - propinquity explains a minority of endogamy;
  - education sorting is strong.

  What needs care is the interpretation of coefficients as tastes, of sex
  differences, of status vs homophily, and above all of Paper 3's
  coefficient component as a change in assortativity.
