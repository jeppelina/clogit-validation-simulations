# =============================================================================
# Study C: does a decomposition across two cohorts separate a change in
# assortativity from a change in composition?
# Run from the repository root:  Rscript C_change_over_time/run.R
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 3000; N <- 1000
OUT <- "C_change_over_time/results"

SHARES <- list(c(0.9, 0.1), c(0.7, 0.3))               # the minority grows
SCENARIOS <- list(constant = c(1, 1), harden = c(1, 1.5), soften = c(1, 0.5))
MODELS <- list(well_specified = c("same_group", "alt_minority", "edu_diff"),
               without_constant = c("same_group", "edu_diff"))
tastes <- function(b_group) c(same_group = b_group, edu_diff = -1)

permutations <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(permutations(v[-i]), function(p) c(v[i], p))))

# Split P(d2, b2) - P(d2, b1) across terms (Shapley: average marginal
# contribution over all orders in which terms switch from b1 to b2).
shapley <- function(cs, b1, b2) {
  value <- function(S) { b <- b1; b[S] <- b2[S]; predicted_share(cs, b) }
  out <- setNames(numeric(length(b1)), names(b1)); perms <- permutations(names(b1))
  for (p in perms) { S <- character(0); prev <- value(S)
    for (term in p) { S <- c(S, term); cur <- value(S); out[term] <- out[term] + cur - prev; prev <- cur } }
  out / length(perms)
}

market <- function(cohort, b_group, setting, u0, seed) {
  set.seed(seed)
  pop <- make_population(N, SHARES[[cohort]])
  un <- form_unions(pop, make_utilities(pop, tastes(b_group), tastes(b_group), u0), setting)
  s <- un$sides$men
  list(cs = build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100),
       endogamy = mean(s$edf$group[s$ego] == s$adf$group[s$alt]))
}

run_rep <- function(rep, sc, setting, u0s) {
  seeds <- sample.int(.Machine$integer.max, 2)
  m1 <- market(1, sc[1], setting, u0s[["c1"]], seeds[1])
  m2 <- market(2, sc[2], setting, u0s[["c2"]], seeds[2])
  m2cf <- market(2, sc[1], setting, u0s[["c2cf"]], seeds[2])   # cohort 2's people, cohort 1's tastes
  rows <- list()
  for (m in names(MODELS)) {
    f1 <- fit_clogit(m1$cs, MODELS[[m]]); f2 <- fit_clogit(m2$cs, MODELS[[m]])
    b1 <- setNames(f1$est, f1$term); b2 <- setNames(f2$est, f2$term)
    p11 <- predicted_share(m1$cs, b1); p21 <- predicted_share(m2$cs, b1); p22 <- predicted_share(m2$cs, b2)
    sh <- shapley(m2$cs, b1, b2)
    rows[[m]] <- data.frame(rep = rep, model = m, observed_change = m2$endogamy - m1$endogamy,
      composition = p21 - p11, coefficients = p22 - p21, part_same_group = sh[["same_group"]],
      part_constant = if ("alt_minority" %in% names(sh)) sh[["alt_minority"]] else NA,
      part_education = sh[["edu_diff"]], true_assortativity = m2$endogamy - m2cf$endogamy,
      b_group_1 = b1[["same_group"]], b_group_2 = b2[["same_group"]])
  }
  do.call(rbind, rows)
}

calibrate <- function(cohort, b_group, setting) setting_u0(setting, function(u)
  make_utilities(make_population(N, SHARES[[cohort]]), tastes(b_group), tastes(b_group), u))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- start_cluster()
clusterExport(clu, c("N", "SHARES", "MODELS", "tastes", "permutations", "shapley", "market", "run_rep"))
res <- list(); idx <- 0
for (scn in names(SCENARIOS)) for (setting in SETTINGS) {
  idx <- idx + 1; sc <- SCENARIOS[[scn]]
  u0s <- c(c1 = calibrate(1, sc[1], setting), c2 = calibrate(2, sc[2], setting), c2cf = calibrate(2, sc[1], setting))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, sc = sc, setting = setting, u0s = u0s))
  r$scenario <- scn; r$setting <- setting; r$seed <- SEED + idx
  res[[idx]] <- r
  write.csv(do.call(rbind, res), file.path(OUT, "decomposition.csv"), row.names = FALSE)
  cat(sprintf("done %s | %s\n", scn, setting))
}
stopCluster(clu)
