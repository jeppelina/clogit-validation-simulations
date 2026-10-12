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

**1. With two-sided matching, the model estimates a property of the pair,
scaled down by search frictions** (Module A).
- With balanced groups the estimates track the sum of both sides'
  coefficients (the pair surplus).
- The scale depends on how many potential partners people actually meet:
  0.96 when everyone meets everyone, about 0.6 when people meet about 17
  before pairing.
- The scale is shared closely, though not exactly, across coefficients.
- **Absolute coefficient sizes say as much about the market as about
  tastes.** Comparisons of relative size are safer.

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

**4. The same-group coefficient depends on group size, so a shrinking
majority looks like hardening** (Modules A and D).
- With identical tastes, larger groups get smaller same-group coefficients
  and smaller groups larger ones. At 90/8/2: majority 0.49, 8% group 1.59,
  relative to the pair surplus.
- When the majority shrinks across cohorts, Paper 3's composition/coefficient
  decomposition reports a coefficient change toward endogamy (+3 to +4
  points of majority endogamy) even when tastes do not change. That is 23
  to 51% of the composition component, depending on taste strength.
- The effect is robust to:
  - which sex proposes;
  - the outside option;
  - market size;
  - which sex is the ego;
  - Paper 3's own alter-category specification.
- Real changes in tastes keep their sign when they are larger than the
  bias, but a real softening is understated by about half.
- **A margin-free measure, the log odds ratio of the couple table, stays
  flat when tastes are flat and moves when they move.** It is a check that
  can be run on the register data.

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
