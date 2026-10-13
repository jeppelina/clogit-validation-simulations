# =============================================================================
# Module A diagnostics for the two failures in the production run
# =============================================================================
# Run from the repo root:  Rscript modules/A_recovery/diag_A.R
# Q1 Group size: under Gale-Shapley the pooled same-group rho falls as the
#    majority grows (0.96 at 70/20/10, 0.86 at 90/8/2). Small-numbers effect
#    (the 2% group is 20 people per sex at N = 1,000) or structural? Which
#    group's coefficient moves? -> N = 1,000 vs 3,000; same-group split by
#    ego group.
# Q2 Segregation: with distance in the model, same-group rho rises with
#    segregation under two-sided matching. Functional form? -> linear
#    ln(distance) vs a natural spline in ln(distance) (5 df).
# =============================================================================

library(parallel); library(splines)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 40; SEED <- 20261030; OUT <- "modules/A_recovery/results"
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
B2 <- c(same_group = 1, edu_diff = -1); B3 <- c(B2, ln_dist = -0.5)

cells <- list(
  list(q = "Q1", label = "70/20/10 N=1000", N = 1000, shares = c(.7, .2, .1), seg = NULL, beta = B2, setting = "gs"),
  list(q = "Q1", label = "90/8/2 N=1000",   N = 1000, shares = c(.9, .08, .02), seg = NULL, beta = B2, setting = "gs"),
  list(q = "Q1", label = "90/8/2 N=3000",   N = 3000, shares = c(.9, .08, .02), seg = NULL, beta = B2, setting = "gs"),
  list(q = "Q1", label = "90/8/2 N=1000",   N = 1000, shares = c(.9, .08, .02), seg = NULL, beta = B2, setting = "seq"),
  list(q = "Q1", label = "90/8/2 N=3000",   N = 3000, shares = c(.9, .08, .02), seg = NULL, beta = B2, setting = "seq"),
  list(q = "Q2", label = "segregation 0.9", N = 1000, shares = c(.7, .2, .1), seg = 0.9, beta = B3, setting = "gs"),
  list(q = "Q2", label = "segregation 0.9", N = 1000, shares = c(.7, .2, .1), seg = 0.9, beta = B3, setting = "seq"),
  list(q = "Q2", label = "segregation 0",   N = 1000, shares = c(.7, .2, .1), seg = 0,   beta = B3, setting = "seq")
)

fit_terms <- function(cs, rhs) {
  fit <- clogit(as.formula(paste("chosen ~", rhs, "+ strata(chooser)")), data = cs, method = "breslow")
  data.frame(term = names(coef(fit)), est = unname(coef(fit)))
}

run_rep <- function(rep, ce, u0) {
  pop <- make_population(ce$N, ce$shares, seg = ce$seg)
  ut <- make_utilities(pop, ce$beta, ce$beta, u0)
  un <- form_unions(pop, ut, SETTINGS[[ce$setting]])
  s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
  out <- list()
  if (ce$q == "Q1") {
    g <- s$edf$group[cs$chooser]
    for (k in 1:3) cs[[paste0("sg_g", k)]] <- cs$same_group * (g == k)
    out$pooled <- cbind(model = "pooled", fit_terms(cs, "same_group + edu_diff"))
    out$bygroup <- cbind(model = "by_ego_group", fit_terms(cs, "sg_g1 + sg_g2 + sg_g3 + edu_diff"))
  } else {
    out$linear <- cbind(model = "linear_dist", fit_terms(cs, "same_group + edu_diff + ln_dist"))
    out$spline <- cbind(model = "spline_dist", fit_terms(cs, "same_group + edu_diff + ns(ln_dist, df = 5)"))
  }
  r <- do.call(rbind, out); r$rep <- rep; r$pair_rate <- un$pair_rate
  r
}

clu <- makeCluster(2)
invisible(clusterEvalQ(clu, { library(splines); for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("run_rep", "fit_terms"))
res <- list()
for (i in seq_along(cells)) {
  ce <- cells[[i]]
  u0 <- setting_u0(SETTINGS[[ce$setting]], function(u)
    make_utilities(make_population(ce$N, ce$shares, seg = ce$seg), ce$beta, ce$beta, u))
  clusterSetRNGStream(clu, SEED + i)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, ce = ce, u0 = u0))
  r$q <- ce$q; r$cell <- ce$label; r$setting <- ce$setting; r$commit <- COMMIT
  res[[i]] <- r
  cat(sprintf("done %s %s %s\n", ce$q, ce$label, ce$setting))
}
stopCluster(clu)
res <- do.call(rbind, res)
write.csv(res, file.path(OUT, "diag_A.csv"), row.names = FALSE)
# rho relative to the pair surplus (2b): same-group terms /2, edu_diff / -2,
# ln_dist / -1. Spline terms are reported as raw estimates.
target <- c(same_group = 2, sg_g1 = 2, sg_g2 = 2, sg_g3 = 2, edu_diff = -2, ln_dist = -1)
res$rho <- ifelse(res$term %in% names(target), res$est / target[res$term], NA)
a <- aggregate(cbind(rho, est) ~ q + cell + setting + model + term, res, mean, na.action = na.pass)
a$mcse <- aggregate(rho ~ q + cell + setting + model + term, res, function(x) sd(x) / sqrt(length(x)), na.action = na.pass)$rho
options(width = 160); print(a[!grepl("ns\\(", a$term), ], digits = 3, row.names = FALSE)
