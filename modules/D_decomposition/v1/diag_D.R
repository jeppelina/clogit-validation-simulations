# =============================================================================
# Module D robustness: is the spurious coefficient component in the null
# scenario (constant tastes, growing minorities) robust?
# Variants under two-sided matching, all null scenario, cohort 1 vs 3:
#   gs (baseline), gs_w (women propose), gs_nooutside (u0 = -Inf, everyone
#   pairs), gs_N2000, seq; women egos alongside men egos; and Paper 3's
#   specification: alter-group indicators (alt_g2, alt_g3) instead of one
#   same-group dummy. Oracle as reference. 100 reps per variant.
# Run from the repo root:  Rscript modules/D_decomposition/diag_D.R
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 100; SEED <- 20261070
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
SHARES <- list(c(0.90, 0.07, 0.03), NULL, c(0.70, 0.20, 0.10))
B <- c(same_group = 1, edu_diff = -1)

variants <- list(
  list(name = "oracle",       st = "oracle", N = 1000, outside = TRUE),
  list(name = "gs",           st = "gs",     N = 1000, outside = TRUE),
  list(name = "gs_w",         st = "gs_w",   N = 1000, outside = TRUE),
  list(name = "gs_nooutside", st = "gs",     N = 1000, outside = FALSE),
  list(name = "gs_N2000",     st = "gs",     N = 2000, outside = TRUE),
  list(name = "seq",          st = "seq",    N = 1000, outside = TRUE)
)

decomp <- function(cs1, cs3, terms) {
  f1 <- fit_clogit(cs1, terms); f3 <- fit_clogit(cs3, terms)
  b1 <- setNames(f1$est, f1$term); b3 <- setNames(f3$est, f3$term)
  p11 <- predicted_share(cs1, b1); p31 <- predicted_share(cs3, b1); p33 <- predicted_share(cs3, b3)
  c(composition = p31 - p11, coefficient = p33 - p31)
}

market <- function(c, v, u0) {
  pop <- make_population(v$N, SHARES[[c]])
  ut <- make_utilities(pop, B, B, u0)
  un <- form_unions(pop, ut, SETTINGS[[v$st]])
  out <- list()
  for (sx in c("men", "women")) {
    s <- un$sides[[sx]]
    keep <- s$edf$group[s$ego] == 1
    cs <- build_choice_sets(s$ego[keep], s$alt[keep], s$edf, s$adf, 100)
    out[[sx]] <- list(cs = cs, endog = mean(s$adf$group[s$alt[keep]] == 1))
  }
  out
}

run_rep <- function(rep, v, u0s) {
  m1 <- market(1, v, u0s[1]); m3 <- market(3, v, u0s[3])
  rows <- list()
  for (sx in c("men", "women")) {
    d <- decomp(m1[[sx]]$cs, m3[[sx]]$cs, c("same_group", "edu_diff"))
    rows[[length(rows) + 1]] <- data.frame(rep = rep, ego_sex = sx, spec = "same_group",
      est_composition = d[["composition"]], est_coefficient = d[["coefficient"]],
      observed_change = m3[[sx]]$endog - m1[[sx]]$endog)
    if (!is.null(m1[[sx]]$cs$alt_g2)) {
      d2 <- decomp(m1[[sx]]$cs, m3[[sx]]$cs, c("alt_g2", "alt_g3", "edu_diff"))
      rows[[length(rows) + 1]] <- data.frame(rep = rep, ego_sex = sx, spec = "alter_groups",
        est_composition = d2[["composition"]], est_coefficient = d2[["coefficient"]],
        observed_change = m3[[sx]]$endog - m1[[sx]]$endog)
    }
  }
  do.call(rbind, rows)
}

clu <- makeCluster(2)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("SHARES", "B", "decomp", "market", "run_rep"))
res <- list()
for (i in seq_along(variants)) {
  v <- variants[[i]]
  u0s <- sapply(c(1, 3), function(c) if (!v$outside) -Inf else
    setting_u0(SETTINGS[[v$st]], function(u) make_utilities(make_population(v$N, SHARES[[c]]), B, B, u)))
  u0s <- c(u0s[1], NA, u0s[2])
  clusterSetRNGStream(clu, SEED + i)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, v = v, u0s = u0s))
  r$variant <- v$name; r$commit <- COMMIT
  res[[i]] <- r
  cat(sprintf("done %s: coefficient component (men, same_group) %+.4f\n", v$name,
              mean(r$est_coefficient[r$ego_sex == "men" & r$spec == "same_group"])))
}
stopCluster(clu)
res <- do.call(rbind, res)
write.csv(res, "modules/D_decomposition/results/diag_D.csv", row.names = FALSE)
a <- aggregate(cbind(est_coefficient, est_composition, observed_change) ~ variant + ego_sex + spec, res, mean)
a$coef_mcse <- aggregate(est_coefficient ~ variant + ego_sex + spec, res, function(x) sd(x) / sqrt(length(x)))$est_coefficient
print(a[order(a$spec, a$ego_sex, a$variant), ], digits = 3, row.names = FALSE)
