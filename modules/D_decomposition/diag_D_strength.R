# =============================================================================
# Module D robustness: does the spurious coefficient component depend on how
# strong the same-group taste is? No-change scenario (constant tastes),
# same-group taste 0.5, 1 and 2, gs and oracle, men egos, 100 reps.
# Run from the repo root.
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 100; N <- 1000; SEED <- 20261100
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
SHARES <- list(c(0.90, 0.07, 0.03), NULL, c(0.70, 0.20, 0.10))
one <- function(coh, b, st, u0) {
  pop <- make_population(N, SHARES[[coh]]); ut <- make_utilities(pop, b, b, u0)
  un <- form_unions(pop, ut, st); s <- un$sides$men; keep <- s$edf$group[s$ego] == 1
  list(cs = build_choice_sets(s$ego[keep], s$alt[keep], s$edf, s$adf, 100),
       endog = mean(s$adf$group[s$alt[keep]] == 1))
}
run_rep <- function(rep, b, st, u0s) {
  m1 <- one(1, b, st, u0s[1]); m3 <- one(3, b, st, u0s[2])
  f1 <- fit_clogit(m1$cs); f3 <- fit_clogit(m3$cs)
  b1 <- setNames(f1$est, f1$term); b3 <- setNames(f3$est, f3$term)
  p11 <- predicted_share(m1$cs, b1); p31 <- predicted_share(m3$cs, b1); p33 <- predicted_share(m3$cs, b3)
  data.frame(rep = rep, est_composition = p31 - p11, est_coefficient = p33 - p31,
             observed_change = m3$endog - m1$endog)
}
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "SHARES", "one", "run_rep"))
res <- list(); i <- 0
for (bg in c(0.5, 1, 2)) for (sn in c("oracle", "gs")) {
  i <- i + 1; b <- c(same_group = bg, edu_diff = -1); st <- SETTINGS[[sn]]
  u0s <- sapply(c(1, 3), function(coh) setting_u0(st, function(u) make_utilities(make_population(N, SHARES[[coh]]), b, b, u)))
  clusterSetRNGStream(clu, SEED + i)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, b = b, st = st, u0s = u0s))
  r$taste <- bg; r$setting <- sn; r$commit <- COMMIT; res[[i]] <- r
  cat(sprintf("taste %.1f %s: coefficient %+.4f, composition %+.4f, observed %+.4f\n", bg, sn,
              mean(r$est_coefficient), mean(r$est_composition), mean(r$observed_change)))
}
stopCluster(clu)
write.csv(do.call(rbind, res), "modules/D_decomposition/results/diag_D_strength.csv", row.names = FALSE)
