# =============================================================================
# Module A production run: blocks A0 to A3, 200 reps per cell
# =============================================================================
# Run from the repo root:  Rscript modules/A_recovery/run_A.R
# Settings (see SPEC.md):
#   oracle   one-sided choice, no outside option        target = b
#   da_full  everyone meets everyone, 85% pairing       target = 2b
#   da_seq   sequential search, budget 50, 10 rounds,   target = 2b
#            85% pairing
# Parallel over reps with a PSOCK cluster (separate processes; forked
# workers crash survival::clogit on macOS). Each cell gets its own RNG
# stream, so results do not depend on the number of workers.
# =============================================================================

library(parallel)
source("R/population.R"); source("R/matching.R"); source("R/estimate.R")

REPS <- 200; SEED <- 20261013; TARGET_PAIR <- 0.85; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/A_recovery/results"
B2 <- c(same_group = 1, edu_diff = -1)
B3 <- c(B2, ln_dist = -0.5)

SETTINGS <- list(
  oracle  = list(mech = "oracle"),
  da_full = list(mech = "da", budget = N, rounds = 1),
  da_seq  = list(mech = "da", budget = 50, rounds = 10)
)

# A cell = block, label, population spec, utility coefficients, and the
# estimation designs to run on each simulated market.
mk_cell <- function(block, label, shares = c(0.7, 0.2, 0.1), seg = NULL,
                    beta = B2, designs)
  list(block = block, label = label, shares = shares, seg = seg, beta = beta,
       designs = designs)

D_BASE <- list(list(design = "uniform_J100", sampler = "uniform", J = 100,
                    terms = names(B2), offset = FALSE))
D_SPACE <- list(
  list(design = "with_dist", sampler = "uniform", J = 100, terms = names(B3), offset = FALSE),
  list(design = "no_dist", sampler = "uniform", J = 100, terms = names(B2), offset = FALSE))
D_SAMPLING <- list(
  list(design = "uniform_J30", sampler = "uniform", J = 30, terms = names(B2), offset = FALSE),
  list(design = "strat_offset", sampler = "strat", J = 100, terms = names(B2), offset = TRUE),
  list(design = "strat_no_offset", sampler = "strat", J = 100, terms = names(B2), offset = FALSE))

CELLS <- c(
  list(mk_cell("A0", "baseline 70/20/10", designs = D_BASE)),
  lapply(list(c(1, 1, 1) / 3, c(0.85, 0.10, 0.05), c(0.90, 0.08, 0.02)), function(sh)
    mk_cell("A1", paste0("shares ", paste(round(100 * sh), collapse = "/")),
            shares = sh, designs = D_BASE)),
  lapply(c(0, 0.5, 0.9), function(sg)
    mk_cell("A2", paste0("segregation ", sg), seg = sg, beta = B3, designs = D_SPACE)),
  list(mk_cell("A3", "sampling designs", designs = D_SAMPLING))
)

# ---- one simulated market, all designs ---------------------------------------
run_rep <- function(rep, cell, st, u0) {
  pop <- make_population(N, cell$shares, seg = cell$seg)
  ut <- make_utilities(pop, cell$beta, cell$beta, u0)
  if (st$mech == "oracle") {
    ch_m <- oracle_choice(ut$U_m, ut$out_m); ch_w <- oracle_choice(t(ut$U_w), ut$out_w)
    men <- which(!is.na(ch_m)); wom <- which(!is.na(ch_w))
    sides <- list(men = list(ego = men, alt = ch_m[men], edf = pop$men, adf = pop$women),
                  women = list(ego = wom, alt = ch_w[wom], edf = pop$women, adf = pop$men))
    pair_rate <- length(men) / N; mult <- 1; met <- N
    endog <- mean(pop$women$group[ch_m[men]][pop$men$group[men] == 1] == 1)
  } else {
    sm <- search_market(ut, st$budget, st$rounds); cp <- sm$couples
    sides <- list(men = list(ego = cp$man, alt = cp$woman, edf = pop$men, adf = pop$women),
                  women = list(ego = cp$woman, alt = cp$man, edf = pop$women, adf = pop$men))
    pair_rate <- nrow(cp) / N; mult <- 2; met <- mean(rowSums(sm$met))
    endog <- mean(pop$women$group[cp$woman][pop$men$group[cp$man] == 1] == 1)
  }
  out <- list()
  for (sx in names(sides)) {
    s <- sides[[sx]]
    cs_cache <- list()
    for (d in cell$designs) {
      key <- paste(d$sampler, d$J)
      if (is.null(cs_cache[[key]]))
        cs_cache[[key]] <- if (d$sampler == "uniform")
          build_choice_sets(s$ego, s$alt, s$edf, s$adf, d$J) else
          build_choice_sets_strat(s$ego, s$alt, s$edf, s$adf, d$J)
      f <- tryCatch(fit_clogit(cs_cache[[key]], d$terms, d$offset),
                    error = function(e) data.frame(term = d$terms, est = NA, se = NA))
      f$design <- d$design; f$ego_sex <- sx
      f$target <- mult * cell$beta[f$term]
      out[[length(out) + 1]] <- f
    }
  }
  r <- do.call(rbind, out)
  r$rep <- rep; r$pair_rate <- pair_rate; r$endog_maj <- endog; r$mean_met <- met
  r
}

# ---- main loop -----------------------------------------------------------------
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(cl, { source("R/population.R"); source("R/matching.R"); source("R/estimate.R"); NULL }))
clusterExport(cl, c("N", "run_rep"))
t0 <- Sys.time(); results <- list(); idx <- 0
for (ci in seq_along(CELLS)) {
  cell <- CELLS[[ci]]
  for (sn in names(SETTINGS)) {
    idx <- idx + 1
    st <- SETTINGS[[sn]]
    u0 <- -Inf
    if (st$mech == "da") {
      mk <- function(u) make_utilities(make_population(N, cell$shares, seg = cell$seg),
                                       cell$beta, cell$beta, u)
      u0 <- calibrate_u0(mk, st$budget, st$rounds, TARGET_PAIR)
    }
    clusterSetRNGStream(cl, SEED + idx)
    reps <- parLapply(cl, seq_len(REPS), run_rep, cell = cell, st = st, u0 = u0)
    r <- do.call(rbind, reps)
    r$block <- cell$block; r$cell <- cell$label; r$setting <- sn; r$u0 <- u0
    r$commit <- COMMIT; r$seed <- SEED + idx
    results[[idx]] <- r
    all <- do.call(rbind, results)
    for (b in unique(all$block))
      write.csv(all[all$block == b, ], file.path(OUT, paste0(b, ".csv")), row.names = FALSE)
    cat(sprintf("[%5.1f min] done %s | %s | %s (u0 = %.2f, pairing %.2f)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                cell$block, cell$label, sn, u0, mean(r$pair_rate)))
  }
}
stopCluster(cl)
cat("Done.\n")
