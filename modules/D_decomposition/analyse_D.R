# =============================================================================
# Module D: summary, checks, figure. Run from the repo root after run_D.R.
# =============================================================================
library(ggplot2)
OUT <- "modules/D_decomposition/results"
r <- read.csv(file.path(OUT, "D.csv"))
r$d_bg <- r$bg_hat2 - r$bg_hat1
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
vars <- c("observed_change", "composition", "coefficients", "part_same_group", "part_constant",
          "part_education", "true_assortativity", "d_bg")
s <- do.call(rbind, lapply(split(r, r[, c("scenario", "setting", "model")], drop = TRUE), function(d) {
  x <- data.frame(scenario = d$scenario[1], setting = d$setting[1], model = d$model[1])
  for (v in vars) { x[[v]] <- mean(d[[v]], na.rm = TRUE); x[[paste0(v, "_mcse")]] <- mcse(d[[v]]) }
  x$sign_agree <- if (all(d$true_assortativity == 0)) NA else
    mean(sign(d$part_same_group) == sign(d$true_assortativity))
  x }))
s <- s[order(match(s$scenario, c("constant", "harden", "soften")), match(s$setting, c("oracle", "gs")), s$model), ]
rownames(s) <- NULL
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)
options(width = 200)
print(s[, c("scenario", "setting", "model", "observed_change", "composition", "part_same_group",
            "part_constant", "part_education", "true_assortativity", "d_bg", "sign_agree")], digits = 3, row.names = FALSE)

# checks (SPEC)
w <- s[s$model == "well_specified", ]
c1 <- w[w$scenario == "constant", ]
chk <- rbind(
  data.frame(check = paste("constant", c1$setting, "same-group part = 0"), value = c1$part_same_group,
             bound = pmax(0.005, 2.5 * c1$part_same_group_mcse)),
  data.frame(check = paste("constant", c1$setting, "change in same-group coefficient = 0"), value = c1$d_bg,
             bound = pmax(0.005, 2.5 * c1$d_bg_mcse)))
chk$pass <- abs(chk$value) < chk$bound
c2 <- w[w$scenario != "constant", ]
chk <- rbind(chk, data.frame(check = paste(c2$scenario, c2$setting, "sign agreement >= 95%"),
                             value = c2$sign_agree, bound = 0.95, pass = c2$sign_agree >= 0.95))
write.csv(chk, file.path(OUT, "checks.csv"), row.names = FALSE)
cat("\n== Checks ==\n"); print(chk, digits = 3, row.names = FALSE)

# figure: same-group part (estimated) vs true assortativity effect
long <- rbind(
  data.frame(s[, c("scenario", "setting", "model")], what = "Estimated: same-group part",
             v = s$part_same_group, se = s$part_same_group_mcse),
  data.frame(s[s$model == "well_specified", c("scenario", "setting", "model")], what = "True assortativity effect",
             v = s$true_assortativity[s$model == "well_specified"], se = s$true_assortativity_mcse[s$model == "well_specified"]))
long$what <- ifelse(long$what == "True assortativity effect", "True effect (counterfactual)",
                    ifelse(long$model == "well_specified", "Estimated, well-specified", "Estimated, without constant"))
long <- unique(long[, c("scenario", "setting", "what", "v", "se")])
long$what <- factor(long$what, c("True effect (counterfactual)", "Estimated, well-specified", "Estimated, without constant"))
long$scenario <- factor(long$scenario, c("constant", "harden", "soften"),
                        c("Tastes constant", "Tastes harden (1 to 1.5)", "Tastes soften (1 to 0.5)"))
long$setting <- factor(long$setting, c("oracle", "gs"), c("One-sided", "Two-sided (Gale-Shapley)"))
p <- ggplot(long, aes(scenario, 100 * v, fill = what)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = 100 * (v - 1.96 * se), ymax = 100 * (v + 1.96 * se)), position = position_dodge(0.8), width = 0) +
  facet_wrap(~ setting) + scale_fill_manual(values = c("grey40", "#0072B2", "#E69F00")) +
  labs(x = NULL, y = "Change in endogamy share due to assortativity (points)", fill = NULL,
       title = "Separating assortativity from composition across two cohorts",
       subtitle = "The minority grows from 10% to 30%. 200 reps per cell.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "fig_decomposition.pdf"), p, width = 10, height = 4.2)
cat("\nWrote summary.csv, checks.csv, fig_decomposition.pdf\n")
