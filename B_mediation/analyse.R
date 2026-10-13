# =============================================================================
# Study B: bias = researcher minus truth, on the same simulated markets.
# Run from the repository root after B_mediation/run.R.
# =============================================================================
library(ggplot2)
# Save a figure as PDF and as PNG (the PNG is shown in the README).
save_fig <- function(file, plot, width, height) {
  ggsave(file, plot, width = width, height = height)
  ggsave(sub("\\.pdf$", ".png", file), plot, width = width, height = height, dpi = 150, bg = "white")
}
OUT <- "B_mediation/results"
r <- read.csv(file.path(OUT, "mediation.csv"))
key <- c("b_setting", "gamma", "setting", "rep")
m <- merge(r[r$model == "researcher", ], r[r$model == "truth", ], by = key, suffixes = c("", "_truth"))
m$bias_setting <- m$b_mediator - m$b_mediator_truth
m$bias_homogamy <- m$b_full - m$b_full_truth
m$bias_share <- m$share - m$share_truth
mcse <- function(x) sd(x) / sqrt(length(x))
s <- do.call(rbind, lapply(split(m, m[, c("b_setting", "gamma", "setting")], drop = TRUE), function(d)
  data.frame(b_setting = d$b_setting[1], gamma = d$gamma[1], setting = d$setting[1],
             setting_coef = mean(d$b_mediator), setting_coef_truth = mean(d$b_mediator_truth),
             homogamy_coef = mean(d$b_full), homogamy_coef_truth = mean(d$b_full_truth),
             share = mean(d$share), share_truth = mean(d$share_truth),
             bias_setting = mean(d$bias_setting), bias_setting_mcse = mcse(d$bias_setting),
             bias_homogamy = mean(d$bias_homogamy), bias_homogamy_mcse = mcse(d$bias_homogamy),
             bias_share = mean(d$bias_share), bias_share_mcse = mcse(d$bias_share),
             overlap = mean(d$overlap), pairing = mean(d$pairing))))
s <- s[order(s$b_setting, s$gamma, match(s$setting, c("one_sided", "two_sided"))), ]
rownames(s) <- NULL
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)
options(width = 200)
print(s[, c("b_setting", "gamma", "setting", "setting_coef", "setting_coef_truth", "homogamy_coef",
            "homogamy_coef_truth", "share", "share_truth")], digits = 3, row.names = FALSE)

# Coefficient biases in units of the true homogamy coefficient, so that
# one-sided and two-sided markets (whose coefficients differ in scale) are
# comparable. The share is in absolute terms.
u <- abs(s$homogamy_coef_truth)
long <- rbind(
  data.frame(s[, c("b_setting", "gamma", "setting")], what = "Setting coefficient", v = s$bias_setting / u, se = s$bias_setting_mcse / u),
  data.frame(s[, c("b_setting", "gamma", "setting")], what = "Direct homogamy", v = s$bias_homogamy / u, se = s$bias_homogamy_mcse / u),
  data.frame(s[, c("b_setting", "gamma", "setting")], what = "Mediation share", v = s$bias_share, se = s$bias_share_mcse))
long$what <- factor(long$what, c("Setting coefficient", "Direct homogamy", "Mediation share"))
long$b_setting <- factor(long$b_setting, c(0, 1, 2.5), c("Setting without effect", "Weaker setting effect (1)", "Strong setting effect (2.5)"))
long$setting <- factor(long$setting, c("one_sided", "two_sided"), c("One-sided choice", "Two-sided matching (Gale-Shapley)"))
p <- ggplot(long, aes(gamma, v, colour = setting)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() + geom_point(size = 2) + geom_errorbar(aes(ymin = v - 1.96 * se, ymax = v + 1.96 * se), width = 0) +
  facet_grid(what ~ b_setting, scales = "free_y") + scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Confounding strength (weight of the unobserved trait in partner choice)",
       y = "Bias: coefficients in units of the true homogamy coefficient; share in points", colour = NULL,
       title = "Where unobserved confounding lands in a mediation analysis",
       subtitle = "The trait also sorts people into settings. Truth = same model with the trait controlled. 200 reps per cell.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
save_fig(file.path(OUT, "fig_confounding.pdf"), p, width = 10, height = 7)
