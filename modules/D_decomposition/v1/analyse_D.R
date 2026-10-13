# =============================================================================
# Module D: summarise the decomposition against its counterfactual truth.
# Run from the repo root after run_D.R (and diag_D.R for the robustness table).
# =============================================================================
library(ggplot2)
OUT <- "modules/D_decomposition/results"
r <- read.csv(file.path(OUT, "D.csv"))
mcse <- function(x) sd(x) / sqrt(length(x))
s <- do.call(rbind, lapply(split(r, r[, c("scenario", "setting")], drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], setting = d$setting[1],
             obs_change = mean(d$endog3 - d$endog1),
             est_comp = mean(d$est_composition), true_comp = mean(d$true_composition),
             est_coef = mean(d$est_coefficient), est_coef_mcse = mcse(d$est_coefficient),
             true_coef = mean(d$true_coefficient), true_coef_mcse = mcse(d$true_coefficient),
             spurious_coef = mean(d$est_coefficient - d$true_coefficient),
             spurious_mcse = mcse(d$est_coefficient - d$true_coefficient),
             sign_agree = if (all(d$true_coefficient == 0)) NA else
               mean(sign(d$est_coefficient) == sign(d$true_coefficient)),
             bg_hat1 = mean(d$bg_hat1), bg_hat3 = mean(d$bg_hat3),
             endog1 = mean(d$endog1), endog3 = mean(d$endog3))))
rownames(s) <- NULL
s <- s[order(match(s$scenario, c("null", "harden", "soften", "frictions")),
             match(s$setting, c("oracle", "gs", "seq"))), ]
write.csv(s, file.path(OUT, "D_summary.csv"), row.names = FALSE)
options(width = 200)
print(s, digits = 3, row.names = FALSE)

# Figure: estimated vs true coefficient component per scenario x setting
long <- rbind(data.frame(s[, c("scenario", "setting")], what = "true (counterfactual)", v = s$true_coef, se = s$true_coef_mcse),
              data.frame(s[, c("scenario", "setting")], what = "estimated (decomposition)", v = s$est_coef, se = s$est_coef_mcse))
long$setting <- factor(long$setting, c("oracle", "gs", "seq"),
                       c("One-sided (oracle)", "Gale-Shapley", "Sequential search"))
long$scenario <- factor(long$scenario, c("null", "harden", "soften", "frictions"),
                        c("No change\nin tastes", "Tastes\nharden", "Tastes\nsoften", "Frictions fall,\ntastes constant"))
p <- ggplot(long, aes(scenario, 100 * v, fill = what)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = 100 * (v - 1.96 * se), ymax = 100 * (v + 1.96 * se)),
                position = position_dodge(0.75), width = 0) +
  facet_wrap(~ setting, nrow = 1) +
  labs(x = NULL, y = "Coefficient component (percentage points of majority endogamy)", fill = NULL,
       title = "Composition change alone produces a coefficient component under two-sided matching",
       subtitle = "Majority share falls from 90% to 70% across cohorts. Positive = endogamy higher than composition predicts. 200 reps.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "D_components.pdf"), p, width = 11, height = 4.5)

if (file.exists(file.path(OUT, "diag_D.csv"))) {
  d <- read.csv(file.path(OUT, "diag_D.csv"))
  a <- aggregate(cbind(est_coefficient, est_composition, observed_change) ~ variant + ego_sex + spec, d, mean)
  a$coef_mcse <- aggregate(est_coefficient ~ variant + ego_sex + spec, d, mcse)$est_coefficient
  cat("\n== Robustness of the null-scenario coefficient component (truth = 0) ==\n")
  print(a[order(a$spec, a$ego_sex, match(a$variant, c("oracle", "gs", "gs_w", "gs_nooutside", "gs_N2000", "seq"))), ],
        digits = 3, row.names = FALSE)
  write.csv(a, file.path(OUT, "diag_D_summary.csv"), row.names = FALSE)
}
