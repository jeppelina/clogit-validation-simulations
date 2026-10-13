# =============================================================================
# Study C: summary, checks, figure. Run from the repository root after
# C_change_over_time/run.R. Changes are in the endogamy share of men's unions.
# =============================================================================
library(ggplot2)
# Save a figure as PDF and as PNG (the PNG is shown in the README).
save_fig <- function(file, plot, width, height) {
  ggsave(file, plot, width = width, height = height)
  ggsave(sub("\\.pdf$", ".png", file), plot, width = width, height = height, dpi = 150, bg = "white")
}
OUT <- "C_change_over_time/results"
r <- read.csv(file.path(OUT, "decomposition.csv"))
r$change_b_group <- r$b_group_2 - r$b_group_1
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
vars <- c("observed_change", "composition", "coefficients", "part_same_group", "part_constant",
          "part_education", "true_assortativity", "change_b_group")
s <- do.call(rbind, lapply(split(r, r[, c("scenario", "setting", "model")], drop = TRUE), function(d) {
  x <- data.frame(scenario = d$scenario[1], setting = d$setting[1], model = d$model[1])
  for (v in vars) { x[[v]] <- mean(d[[v]], na.rm = TRUE); x[[paste0(v, "_mcse")]] <- mcse(d[[v]]) }
  x$sign_agree <- if (all(d$true_assortativity == 0)) NA else mean(sign(d$part_same_group) == sign(d$true_assortativity))
  x }))
s <- s[order(match(s$scenario, c("constant", "harden", "soften")), match(s$setting, c("one_sided", "two_sided")), s$model), ]
rownames(s) <- NULL
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)
options(width = 200)
print(s[, c("scenario", "setting", "model", "observed_change", "composition", "part_same_group", "part_constant",
            "true_assortativity", "change_b_group", "sign_agree")], digits = 3, row.names = FALSE)

# Checks for the well-specified model.
w <- s[s$model == "well_specified", ]; k <- w[w$scenario == "constant", ]
chk <- rbind(
  data.frame(check = paste("constant tastes,", k$setting, ": same-group part = 0"), value = k$part_same_group,
             bound = pmax(0.005, 2.5 * k$part_same_group_mcse)),
  data.frame(check = paste("constant tastes,", k$setting, ": change in same-group coefficient = 0"), value = k$change_b_group,
             bound = pmax(0.005, 2.5 * k$change_b_group_mcse)))
chk$pass <- abs(chk$value) < chk$bound
h <- w[w$scenario != "constant", ]
chk <- rbind(chk, data.frame(check = paste(h$scenario, h$setting, ": sign agreement >= 95%"),
                             value = h$sign_agree, bound = 0.95, pass = h$sign_agree >= 0.95))
write.csv(chk, file.path(OUT, "checks.csv"), row.names = FALSE)
cat(sprintf("\nChecks passed: %d of %d\n", sum(chk$pass), nrow(chk)))

est <- data.frame(s[, c("scenario", "setting")], what = ifelse(s$model == "well_specified", "Estimated, well-specified", "Estimated, without constant"),
                  v = s$part_same_group, se = s$part_same_group_mcse)
tru <- data.frame(w[, c("scenario", "setting")], what = "True effect (counterfactual)", v = w$true_assortativity, se = w$true_assortativity_mcse)
long <- rbind(tru, est)
long$what <- factor(long$what, c("True effect (counterfactual)", "Estimated, well-specified", "Estimated, without constant"))
long$scenario <- factor(long$scenario, c("constant", "harden", "soften"),
                        c("Tastes constant", "Tastes harden (1 to 1.5)", "Tastes soften (1 to 0.5)"))
long$setting <- factor(long$setting, c("one_sided", "two_sided"), c("One-sided choice", "Two-sided matching (Gale-Shapley)"))
p <- ggplot(long, aes(scenario, 100 * v, fill = what)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = 100 * (v - 1.96 * se), ymax = 100 * (v + 1.96 * se)), position = position_dodge(0.8), width = 0) +
  facet_wrap(~ setting) + scale_fill_manual(values = c("grey40", "#0072B2", "#E69F00")) +
  labs(x = NULL, y = "Change in endogamy share due to assortativity (points)", fill = NULL,
       title = "Separating assortativity from composition across two cohorts",
       subtitle = "The minority grows from 10% to 30%. 200 reps per cell.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
save_fig(file.path(OUT, "fig_decomposition.pdf"), p, width = 10, height = 4.2)
