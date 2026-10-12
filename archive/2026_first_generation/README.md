# First-generation simulations (March to September 2026), archived

This folder holds the simulation work that came before the October 2026
revision. It is kept for the record: the April 2026 kappa appendix, Paper
edits drafted in modules 7 and 8, and earlier notes all cite it. **It is
superseded by `modules/A_recovery` to `modules/D_decomposition` at the
repository root**, which share one core library (`R/`), one validated
partnering rule and a common set of settings.

Paths inside these scripts point to their pre-move locations and will need
adjusting to run again.

## What is here

| Folder / file | What it did | Superseded by |
|---|---|---|
| `module1_group_size.R`, `module1b_identification.R`, `module2_spatial.R`, `module3_mediation.R`, `run_module4.R`, `run_modules1to3.R`, `utils.R`, `results/` | Modules 1 to 4: one-sided vs threshold-bilateral recovery, homophily vs status, spatial confounding, mediation, matching mechanisms. Source of the April 2026 kappa appendix (`generate_appendix_tables.R`). | A, B, C |
| `module5_boundary_matrix/` | Recovery of a K x K boundary matrix; marked experimental, never validated. Tested a log-linear Paper 2 pipeline that Paper 2 no longer uses. | Not carried forward |
| `module6_temporal_composition/` | Decade-based decomposition tests (phases 1 to 8b), including Choo-Siow and pair-symmetric estimators. | D |
| `module7_identification_theory/` | Identification theory (pair surplus), experiments 1 to 6c, draft kappa and paper edits. | B |
| `module8_estimand_search/` | Paper 4 estimand search (KHB "artifact floor", sex gap), Paper 3 control leakage (P9 to P11). | C, B1b |

## Known problems, found in the October 2026 revision

These are why the work was rebuilt rather than extended.

1. **Single-sex population (modules 1 to 4).** `utils.R` draws one pool in
   which anyone can partner anyone; each union is counted from both ends
   ("pair rate 2.00"). The papers estimate two-sex models.
2. **Threshold "bilateral" matching** (modules 1 to 4, 6). Random encounters
   with a fixed or calibrated acceptance threshold, and no competition. The
   results depend heavily on the threshold.
3. **"Deferred acceptance" was not Gale-Shapley.**
   - `utils.R::deferred_acceptance_matching` runs on the single-sex pool,
     samples a few candidates per round, and caps the rounds.
   - Module 7's version (`abstract_dgp.R`), which module 8 calls as
     `gale_shapley`, has three problems:
     - each proposer sees 8 random candidates per round;
     - rounds are capped;
     - taste shocks are redrawn at every comparison.
     So preferences are not fixed and the outcome is not a stable matching.
     Module 8's "the sex gap appears only under competitive matching" rests
     on it.
4. **April appendix tables.**
   - The Module 4 table is labelled "20% calibrated acceptance" but comes
     from the fixed-threshold Part 1 run. Under that threshold every agent
     in the status-maximisation world accepted everyone (random matching).
   - The Mixed-DGP bilateral homophily coefficient is described as
     "inflated by 78%". It is attenuated by 78% (-0.33 against -1.5).
   - "Additivity" holds for the education coefficient but not for the
     workplace coefficient: +7% bias without the confounder, +51% with it.
   - The "15 to 40% at acceptance 0.2 to 0.4" range is actually +37% to
     -16%.
5. **The module 8 "artifact floor"** (a channel with zero effect credited
   with ~0.19 to 0.29 of mediation, about 2 x cor^2).
   - It was measured on a coarse "both partners in the top quarter"
     indicator, while the tastes were continuous in education.
   - With a correctly specified education term, the floor is zero
     (`modules/C_mediation`, C1).
   - The floor is real, but it is a misspecification effect: a coarse
     homogamy measure leaves similarity out, and any education-sorted
     channel absorbs it. Module C's coarse-education addendum shows this
     directly.
6. **Bilateral matching in module 3 silently dropped "level" and
   "interaction" confounders** from the matching utility (a warning was
   added in May 2026). The core stages used the "homophily" confounder and
   are unaffected.

## Do not publish as is

`module8_estimand_search/findings/P11_paper3_real_data_spec_comparison.md`
quotes unpublished Paper 3 coefficients. Remove it, or get Jesper's
sign-off, before anything in this repository is pushed to the public GitHub
remote.
