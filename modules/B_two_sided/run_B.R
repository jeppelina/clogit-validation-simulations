# =============================================================================
# Module B production run: B1a (direction), B1b (same tastes, different
# pools), B2 (homophily vs status). 200 reps per cell. See SPEC.md.
# Run from the repo root:  Rscript modules/B_two_sided/run_B.R
# =============================================================================

library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))

REPS <- 200; SEED <- 20261020; N <- 1000
N_WORKERS <- max(1, detectCores() - 2)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
OUT <- "modules/B_two_sided/results"

# A cell: population arguments, men's and women's tastes, estimated terms.
cell <- function(block, world, beta_m, beta_w, terms = c("same_group", "edu_diff"),
                 pop_args = list(), settings)
  list(block = block, world = world, beta_m = beta_m, beta_w = beta_w,
       terms = terms, pop_args = pop_args, settings = settings)

B1A_SET <- c("oracle", "gs", "gs_w", "seq")
B_SET <- c("oracle", "gs", "seq")
CELLS <- list(
  # B1a: same pair surplus (2 on same_group), different holder
  cell("B1a", "men_only",   c(same_group = 2, edu_diff = -1), c(same_group = 0, edu_diff = -1), settings = B1A_SET),
  cell("B1a", "women_only", c(same_group = 0, edu_diff = -1), c(same_group = 2, edu_diff = -1), settings = B1A_SET),
  cell("B1a", "both",       c(same_group = 1, edu_diff = -1), c(same_group = 1, edu_diff = -1), settings = B1A_SET),
  # B1b: identical tastes, different pools
  cell("B1b", "symmetric", c(same_group = 1, edu_diff = -1), c(same_group = 1, edu_diff = -1), settings = B_SET),
  cell("B1b", "sex_ratio", c(same_group = 1, edu_diff = -1), c(same_group = 1, edu_diff = -1),
       pop_args = list(shares = c(0.64, 0.24, 0.12), shares_w = c(0.76, 0.16, 0.08)), settings = B_SET),
  cell("B1b", "edu_gap",   c(same_group = 1, edu_diff = -1), c(same_group = 1, edu_diff = -1),
       pop_args = list(edu_shift_w = 0.5), settings = B_SET),
  # B2: homophily vs status
  cell("B2", "homophily", c(same_group = 1, edu_diff = -1), c(same_group = 1, edu_diff = -1),
       terms = c("same_group", "edu_diff", "alt_edu"), settings = B_SET),
  cell("B2", "status",    c(same_group = 1, alt_edu = 1),   c(same_group = 1, alt_edu = 1),
       terms = c("same_group", "edu_diff", "alt_edu"), settings = B_SET),
  cell("B2", "mixed",     c(same_group = 1, edu_diff = -1, alt_edu = 1), c(same_group = 1, edu_diff = -1, alt_edu = 1),
       terms = c("same_group", "edu_diff", "alt_edu"), settings = B_SET)
)

run_rep <- function(rep, cel, st, u0) {
  pop <- do.call(make_population, c(list(N = N), cel$pop_args))
  ut <- make_utilities(pop, cel$beta_m, cel$beta_w, u0)
  un <- form_unions(pop, ut, st)
  f <- fit_both_sides(un, cel$terms)
  own <- function(sx, t) { b <- if (sx == "men") cel$beta_m else cel$beta_w; if (t %in% names(b)) b[[t]] else 0 }
  oth <- function(sx, t) { b <- if (sx == "men") cel$beta_w else cel$beta_m; if (t %in% names(b)) b[[t]] else 0 }
  f$b_own <- mapply(own, f$ego_sex, f$term); f$b_other <- mapply(oth, f$ego_sex, f$term)
  f$rep <- rep; f$pair_rate <- un$pair_rate
  # education sorting among couples: correlation of partners' education
  f$edu_cor <- cor(pop$men$edu[un$couples$man], pop$women$edu[un$couples$woman])
  f$endog_maj <- mean(pop$women$group[un$couples$woman][pop$men$group[un$couples$man] == 1] == 1)
  f
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- makeCluster(N_WORKERS)
invisible(clusterEvalQ(clu, { for (f in c("population", "matching", "estimate", "run_helpers"))
  source(file.path("R", paste0(f, ".R"))); NULL }))
clusterExport(clu, c("N", "run_rep"))
t0 <- Sys.time(); res <- list(); idx <- 0
for (cel in CELLS) for (sn in cel$settings) {
  idx <- idx + 1; st <- SETTINGS[[sn]]
  u0 <- setting_u0(st, function(u) {
    pop <- do.call(make_population, c(list(N = N), cel$pop_args))
    make_utilities(pop, cel$beta_m, cel$beta_w, u)
  })
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, cel = cel, st = st, u0 = u0))
  r$block <- cel$block; r$world <- cel$world; r$setting <- sn; r$u0 <- u0
  r$commit <- COMMIT; r$seed <- SEED + idx
  res[[idx]] <- r
  all <- do.call(rbind, res)
  for (b in unique(all$block))
    write.csv(all[all$block == b, ], file.path(OUT, paste0(b, ".csv")), row.names = FALSE)
  cat(sprintf("[%5.1f min] done %s | %s | %s (u0 = %.2f, pairing %.2f)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              cel$block, cel$world, sn, u0, mean(r$pair_rate)))
}
stopCluster(clu)
cat("Done.\n")
