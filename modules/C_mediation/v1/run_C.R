# =============================================================================
# Module C production run: KHB mediation shares under oracle and two-sided
# matching. 200 reps per cell, men as egos. See SPEC.md.
# Run from the repo root:  Rscript modules/C_mediation/run_C.R
# =============================================================================

library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 20261040; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/C_mediation/results"
INST <- list(K = 10, sort = 5)          # 10% overlap among alternatives, share ~0.24 under oracle

cells <- list(
  list(cell = "C1_placebo",    kind = "inst", trait = FALSE,
       beta = c(same_group = 1, edu_diff = -1, same_inst = 0)),
  list(cell = "C2_channel",    kind = "inst", trait = FALSE,
       beta = c(same_group = 1, edu_diff = -1, same_inst = 2.5)),
  list(cell = "C3_confounded", kind = "inst", trait = TRUE,
       beta = c(same_group = 1, edu_diff = -1, same_inst = 2.5, t_diff = -1)),
  list(cell = "C4_seg0.5",     kind = "dist", seg = 0.5,
       beta = c(same_group = 1, edu_diff = -1, ln_dist = -0.5)),
  list(cell = "C4_seg0.9",     kind = "dist", seg = 0.9,
       beta = c(same_group = 1, edu_diff = -1, ln_dist = -0.5))
)

make_pop <- function(ce) {
  if (ce$kind == "inst")
    add_institutions(make_population(N), K = INST$K, sort = INST$sort, trait = ce$trait)
  else make_population(N, seg = ce$seg)
}

# Researcher's KHB, plus (C3 only) the infeasible KHB that controls for the
# unobserved trait: the share the channel would get if t were measured.
run_rep <- function(rep, ce, st, u0) {
  pop <- make_pop(ce)
  ut <- make_utilities(pop, ce$beta, ce$beta, u0)
  un <- form_unions(pop, ut, st)
  s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
  out <- if (ce$kind == "inst") {
    r <- cbind(spec = "as_researcher", khb_share(cs, "edu_diff", "same_inst", "same_group"))
    if (ce$trait) r <- rbind(r, cbind(spec = "trait_controlled",
                       khb_share(cs, "edu_diff", "same_inst", c("same_group", "t_diff"))))
    r
  } else cbind(spec = "as_researcher", khb_share(cs, "same_group", "ln_dist", "edu_diff"))
  out$rep <- rep; out$pair_rate <- un$pair_rate
  out
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "INST", "make_pop", "run_rep"))
t0 <- Sys.time(); res <- list(); idx <- 0
for (ce in cells) for (sn in c("oracle", "gs", "seq")) {
  idx <- idx + 1; st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u) make_utilities(make_pop(ce), ce$beta, ce$beta, u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, ce = ce, st = st, u0 = u0))
  r$cell <- ce$cell; r$setting <- sn; r$u0 <- u0; r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "C.csv"), row.names = FALSE)
  cat(sprintf("[%5.1f min] done %s | %s (pairing %.2f, share %.3f)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), ce$cell, sn,
              mean(r$pair_rate), mean(r$share[r$spec == "as_researcher"])))
}
stopCluster(clu)
cat("Done.\n")
