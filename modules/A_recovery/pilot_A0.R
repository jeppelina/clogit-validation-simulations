# =============================================================================
# Module A, block A0 pilot: what does each mechanism recover, with no threats?
# =============================================================================
# Run from the repo root:  Rscript modules/A_recovery/pilot_A0.R
# Both sexes have identical coefficients, so:
#   oracle target            = own coefficient (b)
#   deferred acceptance      = pair surplus (2b)
# 20 reps per cell: a pilot to see the pattern, not the production run.
# =============================================================================

source("R/population.R"); source("R/matching.R"); source("R/estimate.R")

BETA <- c(same_group = 1, edu_diff = -1)
J <- 100; REPS <- 20; SEED <- 20261012
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)

cells <- rbind(
  data.frame(mech = "oracle", proposer = "none", k = 0, u0 = -Inf, N = 1000),  # k = 0: no meeting network
  expand.grid(mech = "da", proposer = "men", k = c(10, 50, 1000),
              u0 = c(-Inf, 2), N = 1000, stringsAsFactors = FALSE),
  data.frame(mech = "da", proposer = "women", k = 1000, u0 = -Inf, N = 1000),
  data.frame(mech = "da", proposer = "men", k = c(500, 2000), u0 = -Inf,
             N = c(500, 2000))          # complete meetings at other market sizes
)

fit_side <- function(ego_ids, partner_ids, ego_df, alt_df) {
  cs <- build_choice_sets(ego_ids, partner_ids, ego_df, alt_df, J)
  fit_clogit(cs)
}

one_run <- function(cell) {
  pop <- make_population(cell$N)
  ut <- make_utilities(pop, BETA, BETA, cell$u0)
  if (cell$mech == "oracle") {
    ch_m <- oracle_choice(ut$U_m, ut$out_m)            # men choose women
    ch_w <- oracle_choice(t(ut$U_w), ut$out_w)         # women choose men
    men <- which(!is.na(ch_m)); wom <- which(!is.na(ch_w))
    res <- rbind(cbind(ego_sex = "men", fit_side(men, ch_m[men], pop$men, pop$women)),
                 cbind(ego_sex = "women", fit_side(wom, ch_w[wom], pop$women, pop$men)))
    res$target <- BETA[res$term]
    res$pair_rate <- mean(!is.na(ch_m))
    res$endog_maj <- mean(pop$women$group[ch_m[men]] == pop$men$group[men] &
                          pop$men$group[men] == 1) / mean(pop$men$group[men] == 1)
  } else {
    cp <- match_market(ut, cell$k, cell$proposer)
    res <- rbind(cbind(ego_sex = "men", fit_side(cp$man, cp$woman, pop$men, pop$women)),
                 cbind(ego_sex = "women", fit_side(cp$woman, cp$man, pop$women, pop$men)))
    res$target <- 2 * BETA[res$term]
    res$pair_rate <- nrow(cp) / cell$N
    maj <- pop$men$group[cp$man] == 1
    res$endog_maj <- mean(pop$women$group[cp$woman][maj] == 1)
  }
  res
}

set.seed(SEED)
out <- list()
t0 <- Sys.time()
for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  for (r in seq_len(REPS)) {
    res <- one_run(cell)
    out[[length(out) + 1]] <- cbind(cell[rep(1, nrow(res)), ], rep = r, res)
  }
  cat(sprintf("cell %d/%d done (%s, k=%s, u0=%s, N=%d) %.1f min\n", i, nrow(cells),
              cell$mech, cell$k, cell$u0, cell$N,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
raw <- do.call(rbind, out)
raw$commit <- COMMIT; raw$seed <- SEED
dir.create("modules/A_recovery/results", showWarnings = FALSE, recursive = TRUE)
write.csv(raw, "modules/A_recovery/results/pilot_A0_raw.csv", row.names = FALSE)

# Summary: ratio to target, coverage of the target by the 95% CI.
# (Grouping columns must not contain NA: aggregate() drops those rows silently.)
raw$ratio <- raw$est / raw$target
raw$covered <- abs(raw$est - raw$target) <= 1.96 * raw$se
key <- c("mech", "proposer", "k", "u0", "N", "ego_sex", "term")
summ <- aggregate(cbind(est, target, ratio, covered, pair_rate, endog_maj) ~ .,
                  data = raw[, c(key, "est", "target", "ratio", "covered",
                                 "pair_rate", "endog_maj")], FUN = mean,
                  na.action = na.pass)
summ$ratio_sd <- aggregate(ratio ~ ., data = raw[, c(key, "ratio")], FUN = sd,
                           na.action = na.pass)$ratio
write.csv(summ, "modules/A_recovery/results/pilot_A0_summary.csv", row.names = FALSE)
print(summ[order(summ$term, summ$mech, summ$N, summ$k, summ$u0, summ$ego_sex), ],
      digits = 3, row.names = FALSE)
