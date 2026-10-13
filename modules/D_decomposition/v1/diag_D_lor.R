# =============================================================================
# Module D follow-up: is a margin-free measure stable where the clogit
# coefficient is not?
# For each cohort market (majority share 90%, 80%, 70%), compute
#   (a) the majority men's clogit same-group coefficient, and
#   (b) the log odds ratio of the couple table, majority vs minority:
#       log( n(M,M) n(m,m) / ( n(M,m) n(m,M) ) ),
#       which is the Choo-Siow cross-difference (singles terms cancel).
# Scenarios: null (constant tastes) and harden (1, 1.25, 1.5).
# 200 reps per cohort x setting x scenario. Run from the repo root.
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; N <- 1000; SEED <- 20261090
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
SHARES <- list(c(0.90, 0.07, 0.03), c(0.80, 0.13, 0.07), c(0.70, 0.20, 0.10))
SCEN <- list(null = c(1, 1, 1), harden = c(1, 1.25, 1.5))

run_rep <- function(rep, coh, bg, st, u0) {
  b <- c(same_group = bg, edu_diff = -1)
  pop <- make_population(N, SHARES[[coh]])
  ut <- make_utilities(pop, b, b, u0)
  un <- form_unions(pop, ut, st)
  cp <- un$couples
  maj_m <- pop$men$group[cp$man] == 1; maj_w <- pop$women$group[cp$woman] == 1
  tab <- table(factor(maj_m, c(TRUE, FALSE)), factor(maj_w, c(TRUE, FALSE))) + 0.5
  lor <- log(tab[1, 1] * tab[2, 2] / (tab[1, 2] * tab[2, 1]))
  s <- un$sides$men; keep <- s$edf$group[s$ego] == 1
  f <- fit_clogit(build_choice_sets(s$ego[keep], s$alt[keep], s$edf, s$adf, 100))
  data.frame(rep = rep, cohort = coh, lor = lor, bg_hat_majority = f$est[f$term == "same_group"],
             endog_majority = mean(maj_w[maj_m]))
}

clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "SHARES", "run_rep"))
res <- list(); idx <- 0
for (scn in names(SCEN)) for (sn in c("oracle", "gs", "seq")) for (c in 1:3) {
  idx <- idx + 1; st <- SETTINGS[[sn]]; bg <- SCEN[[scn]][c]
  b <- c(same_group = bg, edu_diff = -1)
  u0 <- setting_u0(st, function(u) make_utilities(make_population(N, SHARES[[c]]), b, b, u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, coh = c, bg = bg, st = st, u0 = u0))
  r$scenario <- scn; r$setting <- sn; r$commit <- COMMIT
  res[[idx]] <- r
  cat("done", scn, sn, "cohort", c, "\n")
}
stopCluster(clu)
res <- do.call(rbind, res)
write.csv(res, "modules/D_decomposition/results/diag_D_lor.csv", row.names = FALSE)
a <- aggregate(cbind(lor, bg_hat_majority, endog_majority) ~ scenario + setting + cohort, res, mean)
a$lor_mcse <- aggregate(lor ~ scenario + setting + cohort, res, function(x) sd(x) / sqrt(length(x)))$lor
a$bg_mcse <- aggregate(bg_hat_majority ~ scenario + setting + cohort, res, function(x) sd(x) / sqrt(length(x)))$bg_hat_majority
print(a[order(a$scenario, a$setting, a$cohort), ], digits = 3, row.names = FALSE)
