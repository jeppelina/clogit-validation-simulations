# =============================================================================
# Study B: where does unobserved confounding land in a mediation analysis?
# Run from the repository root:  Rscript B_mediation/run.R
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 2000; N <- 1000
OUT <- "B_mediation/results"

make_pop <- function() add_settings(make_population(N, shares = 1), K = 10, sort = 5)
tastes <- function(b_setting, gamma) c(edu_diff = -1, same_setting = b_setting, t_diff = -gamma)

# The researcher's KHB decomposition, and the same decomposition with the
# unobserved trait controlled (the truth, available only in a simulation).
run_rep <- function(rep, b_setting, gamma, setting, u0) {
  pop <- make_pop(); b <- tastes(b_setting, gamma)
  un <- form_unions(pop, make_utilities(pop, b, b, u0), setting); s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
  r <- rbind(cbind(model = "researcher", khb_share(cs, "edu_diff", "same_setting")),
             cbind(model = "truth", khb_share(cs, "edu_diff", "same_setting", "t_diff")))
  r$rep <- rep; r$pairing <- un$pairing; r$overlap <- mean(cs$same_setting[cs$chosen == 0])
  r
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- start_cluster(); clusterExport(clu, c("N", "make_pop", "tastes", "run_rep"))
res <- list(); idx <- 0
for (b_setting in c(0, 1, 2.5)) for (gamma in c(0, 0.5, 1)) for (setting in SETTINGS) {
  idx <- idx + 1
  u0 <- setting_u0(setting, function(u) make_utilities(make_pop(), tastes(b_setting, gamma), tastes(b_setting, gamma), u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, b_setting = b_setting, gamma = gamma, setting = setting, u0 = u0))
  r$b_setting <- b_setting; r$gamma <- gamma; r$setting <- setting; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "mediation.csv"), row.names = FALSE)
  cat(sprintf("done b_setting = %.1f, gamma = %.1f, %s\n", b_setting, gamma, setting))
}
stopCluster(clu)
