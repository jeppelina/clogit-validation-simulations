# =============================================================================
# Module A (v2) diagnostic: why do same-group terms rise with segregation in
# the full specification (distance and alter-group constants in the model)?
# Hypothesis: local margins. Under segregation each person faces a local
# market whose group mix differs from the national one; global alter-group
# constants cannot absorb that. Test: alter-group constants by the ego's
# region (four quadrants of the unit square).
# Segregation 0 and 0.9, gs and seq (and oracle), 100 reps. Run from repo root.
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 100; N <- 1000; SEED <- 20261120
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
B3 <- c(same_group = 1, edu_diff = -1, ln_dist = -0.5)
DIAG <- paste0("diag_g", 1:3)

run_rep <- function(rep, seg, st, u0) {
  pop <- make_population(N, seg = seg); ut <- make_utilities(pop, B3, B3, u0)
  un <- form_unions(pop, ut, st); s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
  g <- s$edf$group[cs$chooser]
  for (k in 1:3) cs[[paste0("diag_g", k)]] <- cs$same_group * (g == k)
  region <- 1 + (s$edf$x[cs$chooser] > 0.5) + 2 * (s$edf$y[cs$chooser] > 0.5)
  local <- character(0)
  for (k in 2:3) for (q in 2:4) {
    v <- paste0("alt_g", k, "_r", q); cs[[v]] <- cs[[paste0("alt_g", k)]] * (region == q); local <- c(local, v)
  }
  base <- c(DIAG, "alt_g2", "alt_g3", "edu_diff", "ln_dist")
  f1 <- fit_clogit(cs, base); f2 <- fit_clogit(cs, c(base, local))
  rbind(cbind(spec = "full_global_constants", f1[f1$term %in% c(DIAG, "ln_dist"), ]),
        cbind(spec = "full_local_constants", f2[f2$term %in% c(DIAG, "ln_dist"), ]),
        data.frame(spec = "pair_rate", term = "pair_rate", est = un$pair_rate, se = NA))
}
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "B3", "DIAG", "run_rep"))
res <- list(); i <- 0
for (seg in c(0, 0.9)) for (sn in c("oracle", "gs", "seq")) {
  i <- i + 1; st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u) make_utilities(make_population(N, seg = seg), B3, B3, u))
  clusterSetRNGStream(clu, SEED + i)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, seg = seg, st = st, u0 = u0))
  r$segregation <- seg; r$setting <- sn; r$commit <- COMMIT; res[[i]] <- r
  cat("done", seg, sn, "\n")
}
stopCluster(clu)
res <- do.call(rbind, res)
write.csv(res, "modules/A_recovery/results/diag_segregation.csv", row.names = FALSE)
res <- res[res$spec != "pair_rate", ]
res$rho <- res$est / (ifelse(res$setting == "oracle", 1, 2) * c(diag_g1 = 1, diag_g2 = 1, diag_g3 = 1, ln_dist = -0.5)[res$term])
a <- aggregate(rho ~ setting + segregation + spec + term, res, mean)
a$mcse <- aggregate(rho ~ setting + segregation + spec + term, res, function(x) sd(x) / sqrt(length(x)))$rho
print(a[order(a$setting, a$spec, a$term, a$segregation), ], digits = 3, row.names = FALSE)
