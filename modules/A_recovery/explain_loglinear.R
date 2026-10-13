# =============================================================================
# Module A: does the group-size effect survive in (a) a clogit with
# alter-group constants, pooled over all ego groups, and (b) a log-linear
# model of the couple table? Three groups, identical tastes (same group 1,
# education -1), everyone pairs. 20 reps per cell. Run from the repo root.
#   clogit_no_constants : diag_g1 + diag_g2 + diag_g3 + edu_diff (Module A spec)
#   clogit_constants    : same + alt_g2 + alt_g3 (column margin)
#   loglinear           : Poisson on the 3 x 3 couple table,
#                         count ~ row + col + diag_g1 + diag_g2 + diag_g3
# Diagonal estimates are reported relative to the pair surplus (2 x 1 under
# two-sided matching): for the log-linear model the diagonal parameter is
# the full pair term, so it is compared with 2b too.
# =============================================================================
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
set.seed(20261014); b <- c(same_group = 1, edu_diff = -1)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
res <- NULL
for (sh in list(c(.7, .2, .1), c(.85, .10, .05))) for (sn in c("oracle", "gs")) for (r in 1:20) {
  pop <- make_population(1000, sh); ut <- make_utilities(pop, b, b, -Inf)
  un <- form_unions(pop, ut, SETTINGS[[sn]]); s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100); g <- s$edf$group[cs$chooser]
  for (k in 1:3) cs[[paste0("diag_g", k)]] <- cs$same_group * (g == k)
  f0 <- fit_clogit(cs, c("diag_g1", "diag_g2", "diag_g3", "edu_diff"))
  f1 <- fit_clogit(cs, c("diag_g1", "diag_g2", "diag_g3", "alt_g2", "alt_g3", "edu_diff"))
  cp <- if (sn == "oracle") data.frame(man = s$ego, woman = s$alt) else un$couples
  tab <- as.data.frame(table(row = factor(pop$men$group[cp$man], 1:3),
                             col = factor(pop$women$group[cp$woman], 1:3)))
  for (k in 1:3) tab[[paste0("diag_g", k)]] <- as.integer(tab$row == k & tab$col == k)
  ll <- coef(glm(Freq ~ row + col + diag_g1 + diag_g2 + diag_g3, family = poisson, data = tab))
  res <- rbind(res, data.frame(shares = paste(round(100 * sh), collapse = "/"), setting = sn, rep = r,
    model = c("clogit_no_constants", "clogit_constants", "loglinear"),
    g1 = c(f0$est[1], f1$est[1], ll[["diag_g1"]]),
    g2 = c(f0$est[2], f1$est[2], ll[["diag_g2"]]),
    g3 = c(f0$est[3], f1$est[3], ll[["diag_g3"]])))
}
res$commit <- COMMIT
write.csv(res, "modules/A_recovery/results/explain_loglinear.csv", row.names = FALSE)
print(aggregate(cbind(g1, g2, g3) ~ setting + shares + model, res, mean), digits = 3, row.names = FALSE)
print(aggregate(cbind(g1, g2, g3) ~ setting + shares + model, res, function(x) sd(x) / sqrt(length(x))), digits = 2, row.names = FALSE)
