# =============================================================================
# Module A: validation of the partnering mechanism before the full A runs
# =============================================================================
# Run from the repo root:  Rscript modules/A_recovery/validate_search.R
#
# V1 Null:    all coefficients zero. Estimates should centre on 0 with ~95%
#             coverage under every mechanism (no manufactured sorting).
# V2 Grid:    search {one-shot, 10 rounds} x budget {1, 5, 20, 100% of N}
#             x outside option {off, calibrated to 85% pairing}, N = 1,000.
# V3 Scaling: does rho depend on the budget as a share of N or as a number?
#             Default search (10 rounds, 85% pairing) at N = 500 and 2,000,
#             with budget fixed at 5% of N and fixed at 50 people.
# Oracle reference cell included. Targets: oracle = b, two-sided = 2b.
# =============================================================================

source("R/population.R"); source("R/matching.R"); source("R/estimate.R")

BETA <- c(same_group = 1, edu_diff = -1)
NULL_BETA <- c(same_group = 0, edu_diff = 0)
J <- 100; SEED <- 20261013; TARGET_PAIR <- 0.85
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/A_recovery/results"

cell <- function(block, mech, N, budget, rounds, outside, null = FALSE, reps = 20)
  data.frame(block, mech, N, budget, rounds, outside, null, reps,
             stringsAsFactors = FALSE)

cells <- rbind(
  # V1 null
  cell("V1_null", "oracle", 1000, NA, NA, "off", TRUE, 50),
  cell("V1_null", "da", 1000, 1000, 1, "off", TRUE, 50),
  cell("V1_null", "da", 1000, 50, 10, "cal85", TRUE, 50),
  # oracle reference
  cell("V2_grid", "oracle", 1000, NA, NA, "off"),
  # V2 grid
  do.call(rbind, lapply(c(10, 50, 200), function(b)
    do.call(rbind, lapply(c(1, 10), function(r)
      do.call(rbind, lapply(c("off", "cal85"), function(o)
        cell("V2_grid", "da", 1000, b, r, o))))))),
  cell("V2_grid", "da", 1000, 1000, 1, "off"),
  cell("V2_grid", "da", 1000, 1000, 1, "cal85"),
  # V3 scaling (N = 1,000 with budget 50 is in V2)
  cell("V3_scale", "da", 500, 25, 10, "cal85"),    # 5% of N
  cell("V3_scale", "da", 500, 50, 10, "cal85"),    # fixed 50
  cell("V3_scale", "da", 2000, 100, 10, "cal85"),  # 5% of N
  cell("V3_scale", "da", 2000, 50, 10, "cal85")    # fixed 50
)

fit_side <- function(ego_ids, partner_ids, ego_df, alt_df)
  fit_clogit(build_choice_sets(ego_ids, partner_ids, ego_df, alt_df, J))

maj_endog <- function(pop, man, woman) {
  maj <- pop$men$group[man] == 1
  mean(pop$women$group[woman][maj] == 1)
}

run_cell <- function(cl) {
  beta <- if (cl$null) NULL_BETA else BETA
  u0 <- -Inf
  if (cl$outside == "cal85") {
    mk <- function(u) make_utilities(make_population(cl$N), beta, beta, u)
    u0 <- calibrate_u0(mk, cl$budget, cl$rounds, TARGET_PAIR)
  }
  res <- vector("list", cl$reps)
  for (r in seq_len(cl$reps)) {
    pop <- make_population(cl$N)
    ut <- make_utilities(pop, beta, beta, u0)
    if (cl$mech == "oracle") {
      ch_m <- oracle_choice(ut$U_m, ut$out_m)
      ch_w <- oracle_choice(t(ut$U_w), ut$out_w)
      men <- which(!is.na(ch_m)); wom <- which(!is.na(ch_w))
      f <- rbind(cbind(ego_sex = "men", fit_side(men, ch_m[men], pop$men, pop$women)),
                 cbind(ego_sex = "women", fit_side(wom, ch_w[wom], pop$women, pop$men)))
      f$target <- beta[f$term]
      f$pair_rate <- length(men) / cl$N
      f$endog_maj <- maj_endog(pop, men, ch_m[men])
      f$rounds_used <- NA; f$mean_met <- cl$N
    } else {
      sm <- search_market(ut, cl$budget, cl$rounds)
      cp <- sm$couples
      f <- rbind(cbind(ego_sex = "men", fit_side(cp$man, cp$woman, pop$men, pop$women)),
                 cbind(ego_sex = "women", fit_side(cp$woman, cp$man, pop$women, pop$men)))
      f$target <- 2 * beta[f$term]
      f$pair_rate <- nrow(cp) / cl$N
      f$endog_maj <- maj_endog(pop, cp$man, cp$woman)
      f$rounds_used <- sm$rounds_used
      f$mean_met <- mean(rowSums(sm$met))
    }
    f$rep <- r; f$u0 <- u0
    res[[r]] <- f
  }
  cbind(cl[rep(1, sum(sapply(res, nrow))), ], do.call(rbind, res))
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
t0 <- Sys.time(); all <- list()
for (i in seq_len(nrow(cells))) {
  all[[i]] <- run_cell(cells[i, ])
  cat(sprintf("[%5.1f min] done cell %2d/%d  %s %s N=%d budget=%s rounds=%s outside=%s\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), i, nrow(cells),
              cells$block[i], cells$mech[i], cells$N[i], cells$budget[i],
              cells$rounds[i], cells$outside[i]))
  raw <- do.call(rbind, all); raw$commit <- COMMIT; raw$seed <- SEED
  write.csv(raw, file.path(OUT, "validate_search_raw.csv"), row.names = FALSE)
}

# ---- summary (grouping columns must not contain NA: aggregate drops them) ----
raw$budget[is.na(raw$budget)] <- 0; raw$rounds[is.na(raw$rounds)] <- 0
raw$rounds_used[is.na(raw$rounds_used)] <- 0
raw$ratio <- ifelse(raw$target == 0, NA, raw$est / raw$target)
raw$covered <- abs(raw$est - raw$target) <= 1.96 * raw$se
key <- c("block", "mech", "N", "budget", "rounds", "outside", "ego_sex", "term")
m <- aggregate(cbind(est, covered, pair_rate, endog_maj, rounds_used, mean_met, u0) ~ .,
               data = raw[, c(key, "est", "covered", "pair_rate", "endog_maj",
                              "rounds_used", "mean_met", "u0")], FUN = mean)
rr <- raw[!is.na(raw$ratio), c(key, "ratio")]
rm <- aggregate(ratio ~ ., data = rr, FUN = mean)
rs <- aggregate(ratio ~ ., data = rr, FUN = sd); names(rs)[ncol(rs)] <- "ratio_sd"
summ <- merge(merge(m, rm, all.x = TRUE), rs, all.x = TRUE)
write.csv(summ, file.path(OUT, "validate_search_summary.csv"), row.names = FALSE)
cat("Done.\n")
