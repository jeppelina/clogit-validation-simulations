# =============================================================================
# Module C: mediation under unobserved confounding. See SPEC.md.
# Run from the repo root:  Rscript modules/C_mediation/run_C.R
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 20261150; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/C_mediation/results"

make_pop <- function() add_institutions(make_population(N, shares = 1), K = 10, sort = 5, trait = TRUE)
tastes <- function(bz, gamma) c(edu_diff = -1, same_inst = bz, t_diff = -gamma)

run_rep <- function(rep, bz, gamma, st, u0) {
  pop <- make_pop(); b <- tastes(bz, gamma)
  ut <- make_utilities(pop, b, b, u0)
  un <- form_unions(pop, ut, st); s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
  res <- rbind(cbind(model = "researcher", khb_share(cs, "edu_diff", "same_inst")),
               cbind(model = "truth", khb_share(cs, "edu_diff", "same_inst", "t_diff")))
  res$rep <- rep; res$pair_rate <- un$pair_rate
  res$overlap <- mean(cs$same_inst[cs$chosen == 0])
  res
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "make_pop", "tastes", "run_rep"))
t0 <- Sys.time(); res <- list(); idx <- 0
for (bz in c(0, 1, 2.5)) for (gamma in c(0, 0.5, 1)) for (sn in c("oracle", "gs")) {
  idx <- idx + 1; st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u) make_utilities(make_pop(), tastes(bz, gamma), tastes(bz, gamma), u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, bz = bz, gamma = gamma, st = st, u0 = u0))
  r$bz <- bz; r$gamma <- gamma; r$setting <- sn; r$u0 <- u0; r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "C.csv"), row.names = FALSE)
  rr <- r[r$model == "researcher", ]; tt <- r[r$model == "truth", ]
  cat(sprintf("[%5.1f min] done bz=%.1f gamma=%.1f %s: setting coef %.2f (truth %.2f), share %.3f (truth %.3f)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), bz, gamma, sn,
              mean(rr$b_mediator), mean(tt$b_mediator), mean(rr$share), mean(tt$share)))
}
stopCluster(clu)
cat("Done.\n")
