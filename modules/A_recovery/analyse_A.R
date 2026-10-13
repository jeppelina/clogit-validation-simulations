# =============================================================================
# Module A (version 3): summary, checks and one figure per test.
# Run from the repo root after run_A.R.
# =============================================================================
library(ggplot2)
OUT <- "modules/A_recovery/results"
r <- do.call(rbind, lapply(file.path(OUT, paste0("test", 1:4, ".csv")), read.csv))
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
key <- c("test", "world", "setting", "model", "ego_sex", "term")
s <- do.call(rbind, lapply(split(r, r[, key], drop = TRUE), function(d)
  data.frame(d[1, key], est = mean(d$est, na.rm = TRUE), mcse = mcse(d$est), truth = d$truth[1],
             n_na = sum(is.na(d$est)), pair_rate = mean(d$pair_rate), edu_cor = mean(d$edu_cor))))
rownames(s) <- NULL
s$rho <- ifelse(!is.na(s$truth) & s$truth != 0, s$est / s$truth, NA)
s$rho_mcse <- ifelse(!is.na(s$rho), s$mcse / abs(s$truth), NA)
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)
SET <- c(oracle = "One-sided", gs = "Two-sided (Gale-Shapley)")
bound <- function(m1, m2 = 0) pmax(0.05, 2.5 * sqrt(m1^2 + m2^2))

# ---- checks ----
chk <- list()
t1 <- s[s$test == "test1" & s$model == "well_specified" & s$term == "same_group", ]
for (st in names(SET)) {
  ref <- t1[t1$setting == st & t1$world == "minority 50%", ]
  for (w in c("minority 20%", "minority 5%")) {
    x <- t1[t1$setting == st & t1$world == w, ]
    chk[[length(chk) + 1]] <- data.frame(test = "test1", check = paste(st, w, "vs 50%"),
      value = x$rho - ref$rho, bound = bound(x$rho_mcse, ref$rho_mcse))
  }
}
for (x in split(t1[t1$setting == "oracle", ], t1$world[t1$setting == "oracle"]))
  chk[[length(chk) + 1]] <- data.frame(test = "test1", check = paste("oracle", x$world, "vs 1"),
                                       value = x$rho - 1, bound = bound(x$rho_mcse))
t2 <- s[s$test == "test2" & s$model == "with_distance" & s$term == "same_group", ]
for (i in seq_len(nrow(t2))) {
  x <- t2[i, ]
  if (x$setting == "oracle") chk[[length(chk) + 1]] <- data.frame(test = "test2",
      check = paste("oracle", x$world, "vs 1"), value = x$rho - 1, bound = bound(x$rho_mcse))
  else if (x$world != "mixed") { ref <- t2[t2$setting == "gs" & t2$world == "mixed", ]
    chk[[length(chk) + 1]] <- data.frame(test = "test2", check = paste("gs", x$world, "vs mixed"),
      value = x$rho - ref$rho, bound = bound(x$rho_mcse, ref$rho_mcse)) }
}
t3 <- s[s$test == "test3" & s$term == "same_group", ]
for (w in unique(t3$world)) for (st in names(SET)) {
  m <- t3[t3$world == w & t3$setting == st & t3$ego_sex == "men", ]
  f <- t3[t3$world == w & t3$setting == st & t3$ego_sex == "women", ]
  if (st == "oracle") { for (x in list(m, f)) chk[[length(chk) + 1]] <- data.frame(test = "test3",
      check = paste("oracle", w, x$ego_sex, "vs own taste"), value = x$est - x$truth, bound = bound(x$mcse)) }
  else chk[[length(chk) + 1]] <- data.frame(test = "test3", check = paste("gs", w, "men minus women"),
      value = m$est - f$est, bound = bound(m$mcse, f$mcse))
}
chk <- do.call(rbind, chk); chk$pass <- abs(chk$value) < chk$bound
write.csv(chk, file.path(OUT, "checks.csv"), row.names = FALSE)
options(width = 180)
print(chk, digits = 3, row.names = FALSE)

# ---- tables printed for the report ----
cat("\n== Test 1 ==\n"); print(s[s$test == "test1" & s$term == "same_group", c("world", "setting", "model", "est", "mcse", "rho")], digits = 3, row.names = FALSE)
cat("\n== Test 2 ==\n"); print(s[s$test == "test2" & s$term %in% c("same_group", "ln_dist"), c("world", "setting", "model", "term", "est", "mcse", "rho")], digits = 3, row.names = FALSE)
cat("\n== Test 3 ==\n"); print(s[s$test == "test3" & s$term == "same_group", c("world", "setting", "ego_sex", "est", "mcse", "truth")], digits = 3, row.names = FALSE)
cat("\n== Test 4 ==\n"); print(s[s$test == "test4" & s$term %in% c("edu_diff", "alt_edu"), c("world", "setting", "term", "est", "mcse", "truth", "edu_cor")], digits = 3, row.names = FALSE)

th <- theme_minimal(base_size = 10) + theme(legend.position = "bottom")
# Test 1
d <- s[s$test == "test1" & s$term == "same_group", ]
d$share <- factor(d$world, c("minority 50%", "minority 20%", "minority 5%"), c("50%", "20%", "5%"))
d$model <- factor(d$model, c("well_specified", "no_constant"), c("With minority-alter constant", "Without"))
d$setting <- factor(SET[d$setting], SET)
ggsave(file.path(OUT, "fig_test1_group_size.pdf"), width = 8, height = 3.8,
  ggplot(d, aes(share, rho, colour = model, group = model)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_line() + geom_point(size = 2) + geom_errorbar(aes(ymin = rho - 1.96 * rho_mcse, ymax = rho + 1.96 * rho_mcse), width = 0) +
    facet_wrap(~ setting) + labs(x = "Minority share", y = "Same-group estimate / target", colour = NULL,
    title = "Test 1. Group size", subtitle = "Identical tastes. Target: true taste (one-sided), pair surplus (two-sided).") + th)
# Test 2
d <- s[s$test == "test2" & s$term == "same_group", ]
d$world <- factor(d$world, c("mixed", "sides moderate", "sides strong", "enclaves moderate", "enclaves strong"),
                  c("mixed", "sides\nmoderate", "sides\nstrong", "enclaves\nmoderate", "enclaves\nstrong"))
d$model <- factor(d$model, c("with_distance", "without_distance"), c("With ln(distance)", "Without"))
d$setting <- factor(SET[d$setting], SET)
ggsave(file.path(OUT, "fig_test2_segregation.pdf"), width = 9, height = 3.8,
  ggplot(d, aes(world, rho, colour = model)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_point(size = 2, position = position_dodge(0.4)) +
    geom_errorbar(aes(ymin = rho - 1.96 * rho_mcse, ymax = rho + 1.96 * rho_mcse), width = 0, position = position_dodge(0.4)) +
    facet_wrap(~ setting) + labs(x = NULL, y = "Same-group estimate / target", colour = NULL,
    title = "Test 2. Segregation", subtitle = "Propinquity in tastes (-0.5 x ln distance). Layouts in layouts.pdf.") + th)
# Test 3
d <- s[s$test == "test3" & s$term == "same_group", ]
d$world <- factor(d$world, c("only men care", "both care", "only women care", "sex ratio"),
                  c("only men\ncare (2, 0)", "both care\n(1, 1)", "only women\ncare (0, 2)", "minority\n60% male (1, 1)"))
d$setting <- factor(SET[d$setting], SET)
ggsave(file.path(OUT, "fig_test3_sex.pdf"), width = 9, height = 3.8,
  ggplot(d, aes(world, est, colour = ego_sex)) + geom_point(size = 2.5, position = position_dodge(0.4)) +
    geom_errorbar(aes(ymin = est - 1.96 * mcse, ymax = est + 1.96 * mcse), width = 0, position = position_dodge(0.4)) +
    geom_point(aes(y = truth), shape = 4, size = 3, position = position_dodge(0.4)) +
    facet_wrap(~ setting) + labs(x = NULL, y = "Same-group estimate (x = target)", colour = "Model fitted for",
    title = "Test 3. Sex differences", subtitle = "Who holds the same-group taste (men, women). Target: own taste (one-sided), pair surplus (two-sided).") + th)
# Test 4
d <- s[s$test == "test4" & s$term %in% c("edu_diff", "alt_edu"), ]
d$term <- factor(d$term, c("edu_diff", "alt_edu"), c("abs(edu difference) (homophily)", "partner's education (status)"))
d$world <- factor(d$world, c("homophily", "status", "both"))
d$setting <- factor(SET[d$setting], SET)
ggsave(file.path(OUT, "fig_test4_status.pdf"), width = 9, height = 3.8,
  ggplot(d, aes(world, est, fill = term)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_point(aes(y = truth), shape = 4, size = 3, position = position_dodge(0.75)) +
    facet_wrap(~ setting) + labs(x = "True tastes", y = "Estimate (x = target)", fill = NULL,
    title = "Test 4. Status vs homophily", subtitle = "Homophily: -1 x abs(edu difference). Status: +2 x partner's education.") + th)
cat("\nWrote summary.csv, checks.csv and four figures.\n")
