# =============================================================================
# Module C addendum: the same institution channel, with education entered as
# a coarse "same tercile" indicator instead of the continuous difference the
# tastes actually use. Tests whether an education-sorted channel absorbs the
# similarity a coarse homogamy measure leaves out (Module 8's "floor").
# 200 reps per cell, men as egos. Run from the repo root.
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 20261080; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/C_mediation/results"
INST <- list(K = 10, sort = 5)
cells <- list(
  list(cell = "C1_placebo_coarse", beta = c(same_group = 1, edu_diff = -1, same_inst = 0)),
  list(cell = "C2_channel_coarse", beta = c(same_group = 1, edu_diff = -1, same_inst = 2.5)))
make_pop <- function() add_institutions(make_population(N), K = INST$K, sort = INST$sort)
run_rep <- function(rep, ce, st, u0) {
  pop <- make_pop(); ut <- make_utilities(pop, ce$beta, ce$beta, u0)
  un <- form_unions(pop, ut, st); s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
  out <- rbind(cbind(spec = "coarse_edu", khb_share(cs, "same_edu3", "same_inst", "same_group")),
               cbind(spec = "continuous_edu", khb_share(cs, "edu_diff", "same_inst", "same_group")))
  out$rep <- rep; out$pair_rate <- un$pair_rate
  out
}
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "INST", "make_pop", "run_rep"))
res <- list(); idx <- 0
for (ce in cells) for (sn in c("oracle", "gs", "seq")) {
  idx <- idx + 1; st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u) make_utilities(make_pop(), ce$beta, ce$beta, u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, ce = ce, st = st, u0 = u0))
  r$cell <- ce$cell; r$setting <- sn; r$u0 <- u0; r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "C_coarse.csv"), row.names = FALSE)
  cat(sprintf("done %s | %s: coarse share %.3f, continuous share %.3f\n", ce$cell, sn,
              mean(r$share[r$spec == "coarse_edu"]), mean(r$share[r$spec == "continuous_edu"])))
}
stopCluster(clu)
cat("Done.\n")
