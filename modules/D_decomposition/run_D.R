# =============================================================================
# Module D: change over two cohorts. See SPEC.md.
# Run from the repo root:  Rscript modules/D_decomposition/run_D.R
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 20261160; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/D_decomposition/results"

SHARES <- list(c(0.9, 0.1), c(0.7, 0.3))
SCEN <- list(constant = c(1, 1), harden = c(1, 1.5), soften = c(1, 0.5))
MODELS <- list(well_specified = c("same_group", "alt_g2", "edu_diff"),
               without_constant = c("same_group", "edu_diff"))
tastes <- function(bg) c(same_group = bg, edu_diff = -1)

all_perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(all_perms(v[-i]), function(p) c(v[i], p))))

# Shapley split of P(d2, b2) - P(d2, b1) across terms.
shapley <- function(cs, b1, b2) {
  terms <- names(b1); val <- function(S) { b <- b1; b[S] <- b2[S]; predicted_share(cs, b) }
  perms <- all_perms(terms); out <- setNames(numeric(length(terms)), terms)
  for (p in perms) { S <- character(0); prev <- val(S)
    for (tm in p) { S <- c(S, tm); cur <- val(S); out[tm] <- out[tm] + (cur - prev); prev <- cur } }
  out / length(perms)
}

market <- function(cohort, bg, st, u0, seed) {
  set.seed(seed)
  pop <- make_population(N, SHARES[[cohort]]); b <- tastes(bg)
  ut <- make_utilities(pop, b, b, u0); un <- form_unions(pop, ut, st); s <- un$sides$men
  list(cs = build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100),
       endog = mean(s$edf$group[s$ego] == s$adf$group[s$alt]))
}

run_rep <- function(rep, sc, st, u0s) {
  seeds <- sample.int(.Machine$integer.max, 2)
  m1 <- market(1, sc[1], st, u0s[["c1"]], seeds[1])
  m2 <- market(2, sc[2], st, u0s[["c2"]], seeds[2])
  m2cf <- market(2, sc[1], st, u0s[["c2cf"]], seeds[2])       # cohort 2's people, cohort 1's tastes
  rows <- list()
  for (m in names(MODELS)) {
    f1 <- fit_clogit(m1$cs, MODELS[[m]]); f2 <- fit_clogit(m2$cs, MODELS[[m]])
    b1 <- setNames(f1$est, f1$term); b2 <- setNames(f2$est, f2$term)
    p11 <- predicted_share(m1$cs, b1); p21 <- predicted_share(m2$cs, b1); p22 <- predicted_share(m2$cs, b2)
    sh <- shapley(m2$cs, b1, b2)
    rows[[m]] <- data.frame(rep = rep, model = m,
      observed_change = m2$endog - m1$endog, composition = p21 - p11, coefficients = p22 - p21,
      part_same_group = sh[["same_group"]],
      part_constant = if ("alt_g2" %in% names(sh)) sh[["alt_g2"]] else NA,
      part_education = sh[["edu_diff"]],
      true_assortativity = m2$endog - m2cf$endog,
      bg_hat1 = b1[["same_group"]], bg_hat2 = b2[["same_group"]])
  }
  do.call(rbind, rows)
}

calib <- function(cohort, bg, st) setting_u0(st, function(u)
  make_utilities(make_population(N, SHARES[[cohort]]), tastes(bg), tastes(bg), u))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "SHARES", "MODELS", "tastes", "all_perms", "shapley", "market", "run_rep"))
t0 <- Sys.time(); res <- list(); idx <- 0
for (scn in names(SCEN)) for (sn in c("oracle", "gs")) {
  idx <- idx + 1; sc <- SCEN[[scn]]; st <- SETTINGS[[sn]]
  u0s <- c(c1 = calib(1, sc[1], st), c2 = calib(2, sc[2], st), c2cf = calib(2, sc[1], st))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, sc = sc, st = st, u0s = u0s))
  r$scenario <- scn; r$setting <- sn; r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "D.csv"), row.names = FALSE)
  w <- r[r$model == "well_specified", ]
  cat(sprintf("[%5.1f min] done %s | %s: same-group part %+.4f (true %+.4f), constant part %+.4f, composition %+.4f\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), scn, sn,
              mean(w$part_same_group), mean(w$true_assortativity), mean(w$part_constant), mean(w$composition)))
}
stopCluster(clu)
cat("Done.\n")
