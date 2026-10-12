# =============================================================================
# Module B2 robustness: does a STRONG status taste turn into apparent
# homophily under Gale-Shapley? Status coefficient 1 (main run) vs 3.
# 100 reps per cell, men as egos. Run from the repo root.
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 100; N <- 1000; SEED <- 20261060
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
terms <- c("same_group", "edu_diff", "alt_edu")
run_rep <- function(rep, b, st, u0) {
  pop <- make_population(N); ut <- make_utilities(pop, b, b, u0)
  un <- form_unions(pop, ut, st); s <- un$sides$men
  f <- fit_clogit(build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100), terms)
  f$rep <- rep; f$edu_cor <- cor(pop$men$edu[un$couples$man], pop$women$edu[un$couples$woman])
  f$pair_rate <- un$pair_rate
  # education of the partnered vs the single (status shows up in who stays single)
  f$edu_partnered_minus_single <- mean(pop$men$edu[un$couples$man]) -
    mean(pop$men$edu[setdiff(seq_len(N), un$couples$man)])
  f
}
clu <- makeCluster(2)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "terms", "run_rep"))
res <- list(); i <- 0
for (bs in c(1, 3)) for (sn in c("oracle", "gs")) {
  i <- i + 1; b <- c(same_group = 1, alt_edu = bs); st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u) make_utilities(make_population(N), b, b, u))
  clusterSetRNGStream(clu, SEED + i)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, b = b, st = st, u0 = u0))
  r$b_status <- bs; r$setting <- sn; r$commit <- COMMIT
  res[[i]] <- r; cat("done", bs, sn, "\n")
}
stopCluster(clu)
res <- do.call(rbind, res)
write.csv(res, "modules/B_two_sided/results/diag_B2.csv", row.names = FALSE)
print(aggregate(cbind(est, edu_cor, pair_rate, edu_partnered_minus_single) ~ b_status + setting + term,
                res, mean), digits = 3)
