# =============================================================================
# Module D production run: composition/coefficient decomposition across three
# cohorts, against a counterfactual truth. 200 reps per cell. See SPEC.md.
# Run from the repo root:  Rscript modules/D_decomposition/run_D.R
# =============================================================================

library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 20261050; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/D_decomposition/results"

SHARES <- list(c(0.90, 0.07, 0.03), c(0.80, 0.13, 0.07), c(0.70, 0.20, 0.10))
SCEN <- list(
  null      = list(bg = c(1, 1, 1),       budget = NULL, settings = c("oracle", "gs", "seq")),
  harden    = list(bg = c(1, 1.25, 1.5),  budget = NULL, settings = c("oracle", "gs", "seq")),
  soften    = list(bg = c(1, 0.75, 0.5),  budget = NULL, settings = c("oracle", "gs", "seq")),
  frictions = list(bg = c(1, 1, 1),       budget = c(25, 50, 100), settings = "seq")
)
beta_of <- function(bg) c(same_group = bg, edu_diff = -1)

# Setting for cohort c (frictions scenario changes the search budget).
st_of <- function(sc, sn, c) {
  st <- SETTINGS[[sn]]
  if (!is.null(sc$budget)) st$budget <- sc$budget[c]
  st
}

# One market: population with cohort c's shares, tastes b, setting st.
# Returns majority men's choice sets and their endogamy.
run_market <- function(c, b, st, u0, seed) {
  set.seed(seed)
  pop <- make_population(N, SHARES[[c]])
  ut <- make_utilities(pop, b, b, u0)
  un <- form_unions(pop, ut, st)
  s <- un$sides$men
  keep <- s$edf$group[s$ego] == 1
  cs <- build_choice_sets(s$ego[keep], s$alt[keep], s$edf, s$adf, 100)
  list(cs = cs, endog = mean(s$adf$group[s$alt[keep]] == 1), pair_rate = un$pair_rate)
}

run_rep <- function(rep, sc, sn, u0s) {
  seeds <- sample.int(.Machine$integer.max, 3)
  m1 <- run_market(1, beta_of(sc$bg[1]), st_of(sc, sn, 1), u0s$c1, seeds[1])
  m3 <- run_market(3, beta_of(sc$bg[3]), st_of(sc, sn, 3), u0s$c3, seeds[3])
  # counterfactual: cohort 3's people, shocks and meetings; cohort 1's tastes
  # and search budget
  m3cf <- run_market(3, beta_of(sc$bg[1]), st_of(sc, sn, 1), u0s$c3cf, seeds[3])
  f1 <- fit_clogit(m1$cs); f3 <- fit_clogit(m3$cs)
  b1 <- setNames(f1$est, f1$term); b3 <- setNames(f3$est, f3$term)
  p11 <- predicted_share(m1$cs, b1); p31 <- predicted_share(m3$cs, b1)
  p33 <- predicted_share(m3$cs, b3)
  data.frame(rep = rep,
             endog1 = m1$endog, endog3 = m3$endog, endog3_cf = m3cf$endog,
             est_composition = p31 - p11, est_coefficient = p33 - p31,
             true_composition = m3cf$endog - m1$endog,
             true_coefficient = m3$endog - m3cf$endog,
             bg_hat1 = b1[["same_group"]], bg_hat3 = b3[["same_group"]],
             be_hat1 = b1[["edu_diff"]], be_hat3 = b3[["edu_diff"]],
             pair1 = m1$pair_rate, pair3 = m3$pair_rate)
}

calib <- function(c, b, st) setting_u0(st, function(u)
  make_utilities(make_population(N, SHARES[[c]]), b, b, u))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "SHARES", "beta_of", "st_of", "run_market", "run_rep"))
t0 <- Sys.time(); res <- list(); idx <- 0
for (scn in names(SCEN)) for (sn in SCEN[[scn]]$settings) {
  idx <- idx + 1; sc <- SCEN[[scn]]
  u0s <- list(c1 = calib(1, beta_of(sc$bg[1]), st_of(sc, sn, 1)),
              c3 = calib(3, beta_of(sc$bg[3]), st_of(sc, sn, 3)),
              c3cf = calib(3, beta_of(sc$bg[1]), st_of(sc, sn, 1)))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, sc = sc, sn = sn, u0s = u0s))
  r$scenario <- scn; r$setting <- sn
  r$u0_c1 <- u0s$c1; r$u0_c3 <- u0s$c3; r$u0_c3cf <- u0s$c3cf
  r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "D.csv"), row.names = FALSE)
  cat(sprintf("[%5.1f min] done %s | %s: est coef %+.4f (true %+.4f), est comp %+.4f (true %+.4f)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), scn, sn,
              mean(r$est_coefficient), mean(r$true_coefficient),
              mean(r$est_composition), mean(r$true_composition)))
}
stopCluster(clu)
cat("Done.\n")
