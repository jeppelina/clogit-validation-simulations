# =============================================================================
# Module A (version 2) production run. See SPEC.md.
# Run from the repo root:  Rscript modules/A_recovery/run_A.R
# Every simulated market is fitted with several specifications; men as egos.
# =============================================================================

library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 20261110; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/A_recovery/results"
B2 <- c(same_group = 1, edu_diff = -1)
B3 <- c(B2, ln_dist = -0.5)

cell <- function(block, label, shares = c(0.7, 0.2, 0.1), seg = NULL, beta = B2)
  list(block = block, label = label, shares = shares, seg = seg, beta = beta)
CELLS <- c(
  list(cell("A0", "70/20/10")),
  lapply(list(c(1, 1, 1) / 3, c(.85, .10, .05), c(.90, .08, .02)), function(sh)
    cell("A1", paste(round(100 * sh), collapse = "/"), shares = sh)),
  lapply(c(0, 0.5, 0.9), function(sg) cell("A2", paste("segregation", sg), seg = sg, beta = B3)),
  list(cell("A3", "sampling"))
)

DIAG <- paste0("diag_g", 1:3); CONST <- c("alt_g2", "alt_g3")

# Fit one specification; returns long rows (term, est, se) or NA rows on failure.
fit_spec <- function(cs, terms, spec, design, offset = FALSE) {
  f <- tryCatch(fit_clogit(cs, terms, offset),
                error = function(e) data.frame(term = terms, est = NA_real_, se = NA_real_))
  cbind(spec = spec, design = design, f)
}

run_rep <- function(rep, ce, st, u0) {
  pop <- make_population(N, ce$shares, seg = ce$seg)
  ut <- make_utilities(pop, ce$beta, ce$beta, u0)
  un <- form_unions(pop, ut, st)
  s <- un$sides$men
  dist <- if (!is.null(ce$seg)) "ln_dist" else character(0)
  add_diag <- function(cs) { g <- s$edf$group[cs$chooser]
    for (k in 1:3) cs[[paste0("diag_g", k)]] <- cs$same_group * (g == k); cs }
  out <- list()
  cs <- add_diag(build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100))
  if (ce$block != "A3") {
    out$full <- fit_spec(cs, c(DIAG, CONST, "edu_diff", dist), "full", "uniform_J100")
    out$noc <- fit_spec(cs, c(DIAG, "edu_diff", dist), "no_constants", "uniform_J100")
    maj <- cs[s$edf$group[cs$chooser] == 1, ]
    out$one <- fit_spec(maj, c(CONST, "edu_diff", dist), "one_group", "uniform_J100")
    if (!is.null(ce$seg))
      out$fnd <- fit_spec(cs, c(DIAG, CONST, "edu_diff"), "full_no_distance", "uniform_J100")
    else {
      cp <- un$couples
      tab <- as.data.frame(table(row = factor(pop$men$group[cp$man], 1:3),
                                 col = factor(pop$women$group[cp$woman], 1:3)))
      for (k in 1:3) tab[[paste0("diag_g", k)]] <- as.integer(tab$row == k & tab$col == k)
      ll <- tryCatch(coef(glm(Freq ~ row + col + diag_g1 + diag_g2 + diag_g3, family = poisson, data = tab)),
                     error = function(e) setNames(rep(NA_real_, 3), DIAG))
      out$ll <- data.frame(spec = "loglinear", design = "couple_table", term = DIAG,
                           est = unname(ll[DIAG]), se = NA_real_)
    }
  } else {
    out$full <- fit_spec(cs, c(DIAG, CONST, "edu_diff"), "full", "uniform_J100")
    cs30 <- add_diag(build_choice_sets(s$ego, s$alt, s$edf, s$adf, 30))
    out$j30 <- fit_spec(cs30, c(DIAG, CONST, "edu_diff"), "full", "uniform_J30")
    cst <- add_diag(build_choice_sets_strat(s$ego, s$alt, s$edf, s$adf, 100, 50))
    out$so <- fit_spec(cst, c(DIAG, CONST, "edu_diff"), "full", "strat_offset", offset = TRUE)
    out$sn <- fit_spec(cst, c(DIAG, CONST, "edu_diff"), "full", "strat_no_offset")
  }
  r <- do.call(rbind, out)
  r$rep <- rep; r$pair_rate <- un$pair_rate; r$mean_met <- un$mean_met
  r$n_g3_egos <- sum(s$edf$group[s$ego] == 3)
  r
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "DIAG", "CONST", "fit_spec", "run_rep"))
t0 <- Sys.time(); res <- list(); idx <- 0
for (ce in CELLS) for (sn in c("oracle", "gs", "seq")) {
  idx <- idx + 1; st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u)
    make_utilities(make_population(N, ce$shares, seg = ce$seg), ce$beta, ce$beta, u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, ce = ce, st = st, u0 = u0))
  r$block <- ce$block; r$cell <- ce$label; r$setting <- sn; r$u0 <- u0
  r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  all <- do.call(rbind, res)
  for (b in unique(all$block))
    write.csv(all[all$block == b, ], file.path(OUT, paste0("v2_", b, ".csv")), row.names = FALSE)
  f <- r[r$spec == "full" & r$design == "uniform_J100" & r$term %in% DIAG, ]
  cat(sprintf("[%5.1f min] done %s | %s | %s  full diag means: %s\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), ce$block, ce$label, sn,
              paste(sprintf("%.2f", tapply(f$est, f$term, mean, na.rm = TRUE)), collapse = " / ")))
}
stopCluster(clu)
cat("Done.\n")
