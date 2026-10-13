# =============================================================================
# Module C: where does confounding bias land? Run from the repo root after
# run_C.R. Bias = researcher minus truth (same simulated markets).
# =============================================================================
library(ggplot2)
OUT <- "modules/C_mediation/results"
r <- read.csv(file.path(OUT, "C.csv"))
key <- c("bz", "gamma", "setting", "rep")
rr <- r[r$model == "researcher", ]; tt <- r[r$model == "truth", ]
m <- merge(rr, tt, by = key, suffixes = c("", "_truth"))
m$bias_setting <- m$b_mediator - m$b_mediator_truth
m$bias_education <- m$b_full - m$b_full_truth
m$bias_share <- m$share - m$share_truth
mcse <- function(x) sd(x) / sqrt(length(x))
s <- do.call(rbind, lapply(split(m, m[, c("bz", "gamma", "setting")], drop = TRUE), function(d)
  data.frame(bz = d$bz[1], gamma = d$gamma[1], setting = d$setting[1],
             setting_coef = mean(d$b_mediator), setting_coef_truth = mean(d$b_mediator_truth),
             bias_setting = mean(d$bias_setting), bias_setting_mcse = mcse(d$bias_setting),
             edu_coef = mean(d$b_full), edu_coef_truth = mean(d$b_full_truth),
             bias_education = mean(d$bias_education), bias_education_mcse = mcse(d$bias_education),
             share = mean(d$share), share_truth = mean(d$share_truth),
             bias_share = mean(d$bias_share), bias_share_mcse = mcse(d$bias_share),
             overlap = mean(d$overlap), pair_rate = mean(d$pair_rate))))
s <- s[order(s$bz, s$gamma, match(s$setting, c("oracle", "gs"))), ]
rownames(s) <- NULL
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)
options(width = 200); print(s, digits = 3, row.names = FALSE)

# sanity checks (SPEC)
cat("\n== Sanity: gamma = 0, researcher equals truth (within 2.5 MCSE) ==\n")
z <- s[s$gamma == 0, ]
print(data.frame(z[, c("bz", "setting")],
                 setting_ok = abs(z$bias_setting) < 2.5 * z$bias_setting_mcse + 1e-9,
                 edu_ok = abs(z$bias_education) < 2.5 * z$bias_education_mcse + 1e-9,
                 share_ok = abs(z$bias_share) < 2.5 * z$bias_share_mcse + 1e-9), row.names = FALSE)
cat("\n== Sanity: one-sided truth model recovers b_z and -1 ==\n")
o <- s[s$setting == "oracle", c("bz", "gamma", "setting_coef_truth", "edu_coef_truth")]
print(o, digits = 3, row.names = FALSE)

# Figure: bias in the three quantities. Coefficient biases are expressed in
# units of the true homogamy coefficient (|education coefficient| in the truth
# model), so one-sided and two-sided markets, whose coefficients differ in
# scale by about 2, are comparable. The share is in absolute points.
u <- abs(s$edu_coef_truth)
long <- rbind(
  data.frame(s[, c("bz", "gamma", "setting")], what = "Setting coefficient",
             v = s$bias_setting / u, se = s$bias_setting_mcse / u),
  data.frame(s[, c("bz", "gamma", "setting")], what = "Direct homogamy",
             v = s$bias_education / u, se = s$bias_education_mcse / u),
  data.frame(s[, c("bz", "gamma", "setting")], what = "Mediation share",
             v = s$bias_share, se = s$bias_share_mcse))
long$what <- factor(long$what, c("Setting coefficient", "Direct homogamy", "Mediation share"))
long$bz <- factor(long$bz, c(0, 1, 2.5), c("Placebo setting (b_z = 0)", "Weaker setting (b_z = 1)", "Strong setting (b_z = 2.5)"))
long$setting <- factor(long$setting, c("oracle", "gs"), c("One-sided", "Two-sided (Gale-Shapley)"))
p <- ggplot(long, aes(gamma, v, colour = setting)) +
  geom_hline(yintercept = 0, colour = "grey50") + geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = v - 1.96 * se, ymax = v + 1.96 * se), width = 0) +
  facet_grid(what ~ bz, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Confounding strength (gamma: weight of the unobserved trait in partner choice)",
       y = "Bias: coefficients in units of the true homogamy coefficient; share in points", colour = NULL,
       title = "Where unobserved confounding lands in a mediation analysis",
       subtitle = "The trait also sorts people into settings. Truth = same model with the trait controlled. 200 reps.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "fig_confounding.pdf"), p, width = 10, height = 7)
cat("\nWrote summary.csv and fig_confounding.pdf\n")
