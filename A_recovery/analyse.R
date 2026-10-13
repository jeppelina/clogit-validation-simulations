# =============================================================================
# Study A: summary table, checks and one figure per test.
# Run from the repository root after A_recovery/run.R.
# =============================================================================
library(ggplot2)
# Save a figure as PDF and as PNG (the PNG is shown in the README).
save_fig <- function(file, plot, width, height) {
  ggsave(file, plot, width = width, height = height)
  ggsave(sub("\\.pdf$", ".png", file), plot, width = width, height = height, dpi = 150, bg = "white")
}
OUT <- "A_recovery/results"
r <- do.call(rbind, lapply(file.path(OUT, paste0("test", 1:4, ".csv")), read.csv))
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
key <- c("test", "world", "setting", "model", "ego_sex", "term")
s <- do.call(rbind, lapply(split(r, r[, key], drop = TRUE), function(d)
  data.frame(d[1, key], est = mean(d$est, na.rm = TRUE), mcse = mcse(d$est), truth = d$truth[1],
             n_failed = sum(is.na(d$est)), pairing = mean(d$pairing), edu_cor = mean(d$edu_cor))))
rownames(s) <- NULL
s$ratio <- ifelse(!is.na(s$truth) & s$truth != 0, s$est / s$truth, NA)
s$ratio_mcse <- ifelse(!is.na(s$ratio), s$mcse / abs(s$truth), NA)
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)

# Checks: a deviation passes if it is below max(0.05, 2.5 Monte Carlo SEs).
bound <- function(a, b = 0) pmax(0.05, 2.5 * sqrt(a^2 + b^2))
chk <- list(); add <- function(test, check, value, bnd) chk[[length(chk) + 1]] <<- data.frame(test, check, value, bound = bnd)
t1 <- s[s$test == "test1" & s$model == "well_specified" & s$term == "same_group", ]
for (st in c("one_sided", "two_sided")) {
  ref <- t1[t1$setting == st & t1$world == "minority 50%", ]
  for (w in c("minority 20%", "minority 5%")) { x <- t1[t1$setting == st & t1$world == w, ]
    add("test1", paste(st, w, "vs minority 50%"), x$ratio - ref$ratio, bound(x$ratio_mcse, ref$ratio_mcse)) }
}
for (x in split(t1[t1$setting == "one_sided", ], t1$world[t1$setting == "one_sided"]))
  add("test1", paste("one_sided", x$world, "vs truth"), x$ratio - 1, bound(x$ratio_mcse))
t2 <- s[s$test == "test2" & s$model == "with_distance" & s$term == "same_group", ]
ref2 <- t2[t2$setting == "two_sided" & t2$world == "mixed", ]
for (i in seq_len(nrow(t2))) { x <- t2[i, ]
  if (x$setting == "one_sided") add("test2", paste("one_sided", x$world, "vs truth"), x$ratio - 1, bound(x$ratio_mcse))
  else if (x$world != "mixed") add("test2", paste("two_sided", x$world, "vs mixed"), x$ratio - ref2$ratio, bound(x$ratio_mcse, ref2$ratio_mcse)) }
t3 <- s[s$test == "test3" & s$term == "same_group", ]
for (w in unique(t3$world)) for (st in c("one_sided", "two_sided")) {
  m <- t3[t3$world == w & t3$setting == st & t3$ego_sex == "men", ]; f <- t3[t3$world == w & t3$setting == st & t3$ego_sex == "women", ]
  if (st == "one_sided") for (x in list(m, f)) add("test3", paste("one_sided", w, x$ego_sex, "vs own taste"), x$est - x$truth, bound(x$mcse))
  else add("test3", paste("two_sided", w, "men minus women"), m$est - f$est, bound(m$mcse, f$mcse))
}
chk <- do.call(rbind, chk); chk$pass <- abs(chk$value) < chk$bound
write.csv(chk, file.path(OUT, "checks.csv"), row.names = FALSE)
options(width = 180)
cat(sprintf("Checks passed: %d of %d\n", sum(chk$pass), nrow(chk)))
if (any(!chk$pass)) print(chk[!chk$pass, ], digits = 3, row.names = FALSE)
show <- function(d, cols) print(d[, cols], digits = 3, row.names = FALSE)
cat("\n== Test 1 ==\n"); show(s[s$test == "test1" & s$term == "same_group", ], c("world", "setting", "model", "est", "ratio"))
cat("\n== Test 2 ==\n"); show(s[s$test == "test2" & s$term %in% c("same_group", "ln_dist"), ], c("world", "setting", "model", "term", "est", "ratio"))
cat("\n== Test 3 ==\n"); show(s[s$test == "test3" & s$term == "same_group", ], c("world", "setting", "ego_sex", "est", "truth"))
cat("\n== Test 4 ==\n"); show(s[s$test == "test4" & s$term %in% c("edu_diff", "alt_edu"), ], c("world", "setting", "term", "est", "truth", "edu_cor"))

SET <- c(one_sided = "One-sided choice", two_sided = "Two-sided matching (Gale-Shapley)")
th <- theme_minimal(base_size = 10) + theme(legend.position = "bottom")
eb <- function(y, se, ...) geom_errorbar(aes(ymin = {{ y }} - 1.96 * {{ se }}, ymax = {{ y }} + 1.96 * {{ se }}), width = 0, ...)

d <- s[s$test == "test1" & s$term == "same_group", ]
d$share <- factor(d$world, c("minority 50%", "minority 20%", "minority 5%"), c("50%", "20%", "5%"))
d$model <- factor(d$model, c("well_specified", "no_constant"), c("With minority-alter constant", "Without"))
d$setting <- factor(SET[d$setting], SET)
save_fig(file.path(OUT, "fig_test1_group_size.pdf"), width = 8, height = 3.8, plot =
  ggplot(d, aes(share, ratio, colour = model, group = model)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_line() + geom_point(size = 2) + eb(ratio, ratio_mcse) + facet_wrap(~ setting) +
    labs(x = "Minority share", y = "Same-group estimate / target", colour = NULL, title = "Test 1. Group size",
         subtitle = "Identical tastes. Target: true taste (one-sided), pair surplus (two-sided).") + th)

d <- s[s$test == "test2" & s$term == "same_group", ]
d$world <- factor(d$world, c("mixed", "sides moderate", "sides strong", "enclaves moderate", "enclaves strong"),
                  c("mixed", "sides\nmoderate", "sides\nstrong", "enclaves\nmoderate", "enclaves\nstrong"))
d$model <- factor(d$model, c("with_distance", "without_distance"), c("With ln(distance)", "Without"))
d$setting <- factor(SET[d$setting], SET)
save_fig(file.path(OUT, "fig_test2_segregation.pdf"), width = 9, height = 3.8, plot =
  ggplot(d, aes(world, ratio, colour = model)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_point(size = 2, position = position_dodge(0.4)) + eb(ratio, ratio_mcse, position = position_dodge(0.4)) +
    facet_wrap(~ setting) + labs(x = NULL, y = "Same-group estimate / target", colour = NULL, title = "Test 2. Segregation",
    subtitle = "Propinquity in tastes (-0.5 x ln distance). Layouts in layouts.pdf.") + th)

d <- s[s$test == "test3" & s$term == "same_group", ]
d$world <- factor(d$world, c("only men care", "both care", "only women care", "sex ratio"),
                  c("only men\ncare (2, 0)", "both care\n(1, 1)", "only women\ncare (0, 2)", "minority\n60% male (1, 1)"))
d$setting <- factor(SET[d$setting], SET)
save_fig(file.path(OUT, "fig_test3_sex.pdf"), width = 9, height = 3.8, plot =
  ggplot(d, aes(world, est, colour = ego_sex)) + geom_point(size = 2.5, position = position_dodge(0.4)) +
    eb(est, mcse, position = position_dodge(0.4)) + geom_point(aes(y = truth), shape = 4, size = 3, position = position_dodge(0.4)) +
    facet_wrap(~ setting) + labs(x = NULL, y = "Same-group estimate (x = target)", colour = "Model fitted for",
    title = "Test 3. Sex differences", subtitle = "Same-group taste (men, women). Target: own taste (one-sided), pair surplus (two-sided).") + th)

d <- s[s$test == "test4" & s$term %in% c("edu_diff", "alt_edu"), ]
d$term <- factor(d$term, c("edu_diff", "alt_edu"), c("|edu difference| (homophily)", "partner's education (status)"))
d$world <- factor(d$world, c("homophily", "status", "both"))
d$setting <- factor(SET[d$setting], SET)
save_fig(file.path(OUT, "fig_test4_status.pdf"), width = 9, height = 3.8, plot =
  ggplot(d, aes(world, est, fill = term)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(position = position_dodge(0.75), width = 0.7) + geom_point(aes(y = truth), shape = 4, size = 3, position = position_dodge(0.75)) +
    facet_wrap(~ setting) + labs(x = "True tastes", y = "Estimate (x = target)", fill = NULL, title = "Test 4. Status vs homophily",
    subtitle = "Homophily: -1 x |edu difference|. Status: +2 x partner's education.") + th)
