# =============================================================================
# Module A, block A3 side check: patched fastclogit kernel vs survival::clogit
# =============================================================================
# Run from the repo root:  Rscript modules/A_recovery/check_fastclogit.R
# Needs the patched kernel (fixed 2026-09-11; see KNOWN_ISSUES_fastclogit.md
# in the dissertation root). Set FASTCLOGIT_LIB to its folder; the default is
# Paper 4's copy, which is not part of this public repository.
#
# For each setting, 20 markets; stratified choice sets with the
# McFadden-Manski offset at two strengths of oversampling. Both estimators
# fit the same data. Pass: coefficients agree to 1e-6 and every fastclogit
# fit converges on the primary (gradient) criterion.
# =============================================================================

FASTCLOGIT_LIB <- Sys.getenv("FASTCLOGIT_LIB",
  "../../Paper 4/MONA scripts/lib")
source("R/population.R"); source("R/matching.R"); source("R/estimate.R")
.this_script_dir <- normalizePath(FASTCLOGIT_LIB)
invisible(capture.output(source(file.path(.this_script_dir, "load_fastclogit.R"))))

set.seed(20261014); N <- 1000; REPS <- 20
B <- c(same_group = 1, edu_diff = -1)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
settings <- list(oracle = NULL, da_full = list(budget = N, rounds = 1),
                 da_seq = list(budget = 50, rounds = 10))
rows <- list()
for (sn in names(settings)) for (n_same in c(50, 95)) for (r in seq_len(REPS)) {
  pop <- make_population(N)
  st <- settings[[sn]]
  if (is.null(st)) {
    ut <- make_utilities(pop, B, B, -Inf)
    ch <- oracle_choice(ut$U_m, ut$out_m); ego <- which(!is.na(ch)); alt <- ch[ego]
  } else {
    ut <- make_utilities(pop, B, B, 0)
    cp <- search_market(ut, st$budget, st$rounds)$couples; ego <- cp$man; alt <- cp$woman
  }
  cs <- build_choice_sets_strat(ego, alt, pop$men, pop$women, 100, n_same)
  sv <- fit_clogit(cs, names(B), offset = TRUE)
  fc <- fastclogit(as.matrix(cs[, names(B)]), cs$chosen, cs$chooser,
                   offset = cs$offset_mm)
  rows[[length(rows) + 1]] <- data.frame(
    setting = sn, n_same = n_same, rep = r, term = names(B),
    est_survival = sv$est, est_fastclogit = unname(fc$coefficients[names(B)]),
    se_survival = sv$se, se_fastclogit = unname(fc$se[names(B)]),
    criterion = fc$convergence_criterion, iterations = fc$iterations,
    offset_spread = diff(range(cs$offset_mm)), commit = COMMIT)
}
res <- do.call(rbind, rows)
res$abs_diff <- abs(res$est_survival - res$est_fastclogit)
write.csv(res, "modules/A_recovery/results/A3_fastclogit_check.csv", row.names = FALSE)
cat(sprintf("fits: %d | max |coef diff| = %.2e | max |se diff| = %.2e\n",
            nrow(res) / 2, max(res$abs_diff),
            max(abs(res$se_survival - res$se_fastclogit))))
print(table(res$criterion[res$term == "same_group"], res$setting[res$term == "same_group"]))
cat(sprintf("offset spread within choice sets: %.1f to %.1f\n",
            min(res$offset_spread), max(res$offset_spread)))
