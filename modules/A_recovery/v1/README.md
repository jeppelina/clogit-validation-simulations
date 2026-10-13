# Module A, version 1 (2026-10-13), superseded

The first full Module A run. One specification: a same-group dummy and
education, with no alter-group constants, egos from all groups pooled.

It found the pipeline sound under one-sided choice, but under two-sided
matching the same-group coefficient moved with group size and segregation.
`../explain_group_size.R` and `../explain_loglinear.R` then showed the
reason: without alter-group constants (the column margin), the same-group
term absorbs how attainable each group is in a two-sided market. A logit
with alter-group constants and all ego groups pooled matches a log-linear
model of the couple table and does not move.

Version 2 (`../SPEC.md`) tests the well-specified model directly. Kept here
for the record; scripts reference paths from before the move.
