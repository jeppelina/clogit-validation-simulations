# =============================================================================
# Study A: does a well-specified conditional logit recover what it should?
# Four tests: group size, segregation, sex differences, status vs homophily.
# Run from the repository root:  Rscript A_recovery/run.R
# =============================================================================
library(parallel)
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
REPS <- 200; SEED <- 1000; N <- 1000
OUT <- "A_recovery/results"

WS <- c("same_group", "alt_minority", "edu_diff")       # the well-specified model
tastes <- function(sg = 1, edu = -1, dist = NULL, status = NULL) {
  b <- c(same_group = sg)
  if (!is.null(edu)) b <- c(b, edu_diff = edu)
  if (!is.null(dist)) b <- c(b, ln_dist = dist)
  if (!is.null(status)) b <- c(b, alt_edu = status)
  b
}
cell <- function(test, world, models, pop = list(), bm = tastes(), bw = bm, sexes = "men")
  list(test = test, world = world, models = models, pop = pop, bm = bm, bw = bw, sexes = sexes)

CELLS <- c(
  # Test 1: group size
  lapply(c(0.50, 0.20, 0.05), function(p)
    cell("test1", sprintf("minority %d%%", round(100 * p)),
         list(well_specified = WS, no_constant = c("same_group", "edu_diff")),
         pop = list(shares = c(1 - p, p)))),
  # Test 2: segregation
  lapply(list(c("mixed", 0), c("sides", 0.5), c("sides", 0.9), c("enclaves", 0.5), c("enclaves", 0.9)), function(l)
    cell("test2", if (l[1] == "mixed") "mixed" else paste(l[1], ifelse(l[2] == "0.5", "moderate", "strong")),
         list(with_distance = c(WS, "ln_dist"), without_distance = WS),
         pop = list(layout = l[1], s = as.numeric(l[2])), bm = tastes(dist = -0.5))),
  # Test 3: sex differences
  list(cell("test3", "only men care", list(well_specified = WS), bm = tastes(sg = 2), bw = tastes(sg = 0),
            sexes = c("men", "women")),
       cell("test3", "only women care", list(well_specified = WS), bm = tastes(sg = 0), bw = tastes(sg = 2),
            sexes = c("men", "women")),
       cell("test3", "both care", list(well_specified = WS), sexes = c("men", "women")),
       cell("test3", "sex ratio", list(well_specified = WS),
            pop = list(shares = c(.76, .24), shares_w = c(.84, .16)), sexes = c("men", "women"))),
  # Test 4: status vs homophily
  list(cell("test4", "homophily", list(well_specified = c(WS, "alt_edu"))),
       cell("test4", "status", list(well_specified = c(WS, "alt_edu")), bm = tastes(edu = NULL, status = 2)),
       cell("test4", "both", list(well_specified = c(WS, "alt_edu")), bm = tastes(status = 2)))
)

run_rep <- function(rep, ce, setting, u0) {
  pop <- do.call(make_population, c(list(N = N), ce$pop))
  ut <- make_utilities(pop, ce$bm, ce$bw, u0)
  un <- form_unions(pop, ut, setting)
  out <- list()
  for (sx in ce$sexes) {
    s <- un$sides[[sx]]
    cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100)
    for (m in names(ce$models)) {
      f <- tryCatch(fit_clogit(cs, ce$models[[m]]),
                    error = function(e) data.frame(term = ce$models[[m]], est = NA_real_, se = NA_real_))
      f$model <- m; f$ego_sex <- sx
      out[[length(out) + 1]] <- f
    }
  }
  r <- do.call(rbind, out)
  r$rep <- rep; r$pairing <- un$pairing
  r$edu_cor <- cor(pop$men$edu[un$couples$man], pop$women$edu[un$couples$woman])
  r
}

# Target: the true taste (one-sided) or the pair surplus, the sum of both
# sides' terms (two-sided). For the partner's education, an ego's pair
# surplus is the ego's own status taste (the partner's taste for the ego's
# education is the same across the ego's alternatives). The minority-alter
# constant is 0 one-sided and has no two-sided target.
truth <- function(ce, sx, term, setting) {
  own <- if (sx == "men") ce$bm else ce$bw; oth <- if (sx == "men") ce$bw else ce$bm
  g <- function(b) if (term %in% names(b)) b[[term]] else 0
  if (term == "alt_minority") return(if (setting == "one_sided") 0 else NA)
  if (setting == "one_sided" || term == "alt_edu") g(own) else g(own) + g(oth)
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clu <- start_cluster(); clusterExport(clu, c("N", "run_rep"))
res <- list(); idx <- 0
for (ce in CELLS) for (setting in SETTINGS) {
  idx <- idx + 1
  u0 <- setting_u0(setting, function(u)
    make_utilities(do.call(make_population, c(list(N = N), ce$pop)), ce$bm, ce$bw, u))
  clusterSetRNGStream(clu, SEED + idx)
  r <- do.call(rbind, parLapply(clu, seq_len(REPS), run_rep, ce = ce, setting = setting, u0 = u0))
  r$test <- ce$test; r$world <- ce$world; r$setting <- setting; r$seed <- SEED + idx
  r$truth <- mapply(function(sx, tm) truth(ce, sx, tm, setting), r$ego_sex, r$term)
  res[[idx]] <- r
  all <- do.call(rbind, res)
  for (tt in unique(all$test)) write.csv(all[all$test == tt, ], file.path(OUT, paste0(tt, ".csv")), row.names = FALSE)
  cat(sprintf("done %s | %s | %s\n", ce$test, ce$world, setting))
}
stopCluster(clu)
