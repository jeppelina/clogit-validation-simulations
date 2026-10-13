# =============================================================================
# Module A: why group size moves the same-group coefficient under two-sided
# matching. Two groups, identical tastes, everyone pairs (no outside option).
# Per-group same-group coefficients (majority egos, minority egos) vs the log
# odds ratio of the couple table, at shares 50/50, 70/30, 90/10.
# 20 reps per cell. Run from the repo root.
# =============================================================================
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
set.seed(20261013); b <- c(same_group = 1, edu_diff = -1)
COMMIT <- system("git rev-parse --short HEAD", intern = TRUE)
res <- NULL
for (sh in list(c(.5, .5), c(.7, .3), c(.9, .1))) for (sn in c("oracle", "gs")) for (r in 1:20) {
  pop <- make_population(1000, sh); ut <- make_utilities(pop, b, b, -Inf)
  un <- form_unions(pop, ut, SETTINGS[[sn]]); s <- un$sides$men
  cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, 100); g <- s$edf$group[cs$chooser]
  cs$sg_maj <- cs$same_group * (g == 1); cs$sg_min <- cs$same_group * (g == 2)
  f <- fit_clogit(cs, c("sg_maj", "sg_min", "edu_diff"))
  cp <- if (sn == "oracle") data.frame(man = s$ego, woman = s$alt) else un$couples
  tab <- table(factor(pop$men$group[cp$man], 1:2), factor(pop$women$group[cp$woman], 1:2))
  res <- rbind(res, data.frame(shares = paste(100 * sh, collapse = "/"), setting = sn, rep = r,
    b_majority = f$est[1], b_minority = f$est[2],
    log_odds_ratio = log(tab[1, 1] * tab[2, 2] / (tab[1, 2] * tab[2, 1])),
    endog_majority_men = tab[1, 1] / sum(tab[1, ]), endog_minority_men = tab[2, 2] / sum(tab[2, ]),
    majority_men_with_minority_women = tab[1, 2], minority_men_with_majority_women = tab[2, 1]))
}
res$commit <- COMMIT
write.csv(res, "modules/A_recovery/results/explain_group_size.csv", row.names = FALSE)
a <- aggregate(cbind(b_majority, b_minority, log_odds_ratio, endog_majority_men, endog_minority_men,
                     majority_men_with_minority_women, minority_men_with_majority_women) ~ shares + setting, res, mean)
a$sum_of_coefficients <- a$b_majority + a$b_minority
print(a[order(a$setting, a$shares), ], digits = 3, row.names = FALSE)
