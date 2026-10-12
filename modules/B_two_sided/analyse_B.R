# =============================================================================
# Module B: summarise B1a, B1b, B2. Run from the repo root after run_B.R.
# Writes results/B_summary.csv and results/B1a_direction.pdf
# =============================================================================

library(ggplot2)
OUT <- "modules/B_two_sided/results"
raw <- do.call(rbind, lapply(file.path(OUT, c("B1a.csv", "B1b.csv", "B2.csv")), function(f) {
  d <- read.csv(f); d[, c("block", "world", "setting", "ego_sex", "term", "est", "se",
                          "b_own", "b_other", "rep", "pair_rate", "edu_cor", "endog_maj")]
}))
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
s <- do.call(rbind, lapply(split(raw, raw[, c("block", "world", "setting", "ego_sex", "term")], drop = TRUE),
  function(d) data.frame(d[1, c("block", "world", "setting", "ego_sex", "term")],
                         est = mean(d$est, na.rm = TRUE), mcse = mcse(d$est),
                         b_own = d$b_own[1], b_other = d$b_other[1],
                         n_na = sum(is.na(d$est)), edu_cor = mean(d$edu_cor),
                         pair_rate = mean(d$pair_rate))))
rownames(s) <- NULL
write.csv(s, file.path(OUT, "B_summary.csv"), row.names = FALSE)

wide <- function(d) {
  w <- reshape(d[, c("world", "setting", "term", "ego_sex", "est", "mcse")],
               idvar = c("world", "setting", "term"), timevar = "ego_sex", direction = "wide")
  w$diff_m_minus_w <- w$est.men - w$est.women
  w$diff_mcse <- sqrt(w$mcse.men^2 + w$mcse.women^2)
  w[order(w$term, w$world, w$setting), ]
}
options(width = 180)
cat("\n== B1a: who holds the assortativity (pair surplus on same_group = 2 in every world) ==\n")
print(wide(s[s$block == "B1a", ]), digits = 3, row.names = FALSE)
cat("\n== B1b: identical tastes, different pools ==\n")
print(wide(s[s$block == "B1b", ]), digits = 3, row.names = FALSE)
cat("\n== B2: homophily vs status (men's models; women's are symmetric) ==\n")
b2 <- s[s$block == "B2" & s$ego_sex == "men", c("world", "setting", "term", "est", "mcse", "edu_cor")]
print(reshape(b2[, c("world", "setting", "term", "est")], idvar = c("world", "setting"),
              timevar = "term", direction = "wide"), digits = 3, row.names = FALSE)
cat("\npartners' education correlation by world x setting:\n")
print(aggregate(edu_cor ~ world + setting, s[s$block == "B2", ], mean), digits = 3, row.names = FALSE)

d <- s[s$block == "B1a" & s$term == "same_group", ]
d$world <- factor(d$world, c("men_only", "both", "women_only"),
                  c("Only men care\n(men 2, women 0)", "Both care\n(1 and 1)", "Only women care\n(men 0, women 2)"))
d$setting <- factor(d$setting, c("oracle", "gs", "gs_w", "seq"),
                    c("One-sided (oracle)", "Gale-Shapley, men propose", "Gale-Shapley, women propose",
                      "Sequential search"))
p <- ggplot(d, aes(world, est, colour = ego_sex, group = ego_sex)) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = est - 1.96 * mcse, ymax = est + 1.96 * mcse), width = 0,
                position = position_dodge(0.3)) +
  facet_wrap(~ setting, nrow = 1) +
  labs(x = NULL, y = "Estimated same-group coefficient", colour = "Ego sex",
       title = "Two-sided data cannot say which sex holds the assortativity",
       subtitle = "Pair surplus on same group = 2 in every world. 200 reps per cell.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "B1a_direction.pdf"), p, width = 11, height = 4.2)
cat("\nWrote B_summary.csv and B1a_direction.pdf\n")
