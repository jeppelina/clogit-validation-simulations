# =============================================================================
# Module A: summarise A0 to A3, apply the pass criteria, draw the rho figure
# =============================================================================
# Run from the repo root after run_A.R:  Rscript modules/A_recovery/analyse_A.R
# Writes results/A_summary.csv, results/A_checks.csv, results/A_rho.pdf
#
# rho = estimate / target (target = b for oracle, 2b for DA settings).
# Pass criteria (SPEC.md):
#   oracle: |mean rho - 1| < 0.02 and 95% coverage in [0.93, 0.97]
#   DA:     |mean rho - reference rho| < 0.05, per term, where the reference
#           is A0 for A1 and A3, and segregation 0 for A2 (A2 adds a distance
#           term to utility, so its own no-segregation cell is the like-for-like
#           baseline; the deviation from A0 is reported too).
# Demonstration designs (A2 no_dist, A3 strat_no_offset) are reported but
# not held to the criteria: they are expected to fail.
# =============================================================================

library(ggplot2)
OUT <- "modules/A_recovery/results"
raw <- do.call(rbind, lapply(file.path(OUT, paste0("A", 0:3, ".csv")), read.csv))
raw$rho <- raw$est / raw$target
raw$covered <- abs(raw$est - raw$target) <= 1.96 * raw$se

key <- c("block", "cell", "setting", "design", "term")
s <- do.call(rbind, lapply(split(raw, raw[, key], drop = TRUE), function(d) {
  data.frame(d[1, key], n = sum(!is.na(d$rho)), n_failed = sum(is.na(d$rho)),
             rho = mean(d$rho, na.rm = TRUE), rho_sd = sd(d$rho, na.rm = TRUE),
             rho_mcse = sd(d$rho, na.rm = TRUE) / sqrt(sum(!is.na(d$rho))),
             rho_men = mean(d$rho[d$ego_sex == "men"], na.rm = TRUE),
             rho_women = mean(d$rho[d$ego_sex == "women"], na.rm = TRUE),
             coverage = mean(d$covered, na.rm = TRUE),
             pair_rate = mean(d$pair_rate), endog_maj = mean(d$endog_maj),
             mean_met = mean(d$mean_met), u0 = d$u0[1])
}))
rownames(s) <- NULL

ref_of <- function(block, setting, term, which = c("own", "A0")) {
  which <- match.arg(which)
  r <- if (which == "A0" || block %in% c("A1", "A3"))
    s[s$block == "A0" & s$setting == setting & s$term == term, ] else
    s[s$block == "A2" & s$cell == "segregation 0" & s$design == "with_dist" &
      s$setting == setting & s$term == term, ]
  if (nrow(r)) r$rho[1] else NA
}
s$ref_rho <- mapply(ref_of, s$block, s$setting, s$term, "own")
s$dev_ref <- s$rho - s$ref_rho
s$dev_A0 <- s$rho - mapply(ref_of, s$block, s$setting, s$term, "A0")
s$demonstration <- s$design %in% c("no_dist", "strat_no_offset")
s$pass <- ifelse(s$demonstration, NA,
  ifelse(s$setting == "oracle",
         abs(s$rho - 1) < 0.02 & s$coverage >= 0.93 & s$coverage <= 0.97,
         abs(s$dev_ref) < 0.05 | (is.na(s$dev_ref) & s$block == "A0")))
s <- s[order(s$block, s$cell, s$setting, s$design, s$term), ]
write.csv(s, file.path(OUT, "A_summary.csv"), row.names = FALSE)

# Gap between terms within each cell/setting/design (approximate common scale)
gap <- aggregate(rho ~ block + cell + setting + design,
                 s[s$term %in% c("same_group", "edu_diff"), ],
                 function(x) diff(range(x)))
names(gap)[ncol(gap)] <- "rho_gap_terms"
write.csv(gap, file.path(OUT, "A_term_gap.csv"), row.names = FALSE)

checks <- s[!s$demonstration, c("block", "cell", "setting", "design", "term",
                                "rho", "rho_mcse", "coverage", "ref_rho",
                                "dev_ref", "pass")]
write.csv(checks, file.path(OUT, "A_checks.csv"), row.names = FALSE)

options(width = 200)
cat("\n== All cells ==\n")
print(s[, c("block", "cell", "setting", "design", "term", "rho", "rho_mcse",
            "rho_men", "rho_women", "coverage", "dev_ref", "dev_A0", "pass",
            "pair_rate", "mean_met", "n_failed")], digits = 3, row.names = FALSE)
cat("\n== Failures (non-demonstration) ==\n")
print(checks[!is.na(checks$pass) & !checks$pass, ], digits = 3, row.names = FALSE)
cat("\n== Largest gap between terms ==\n")
print(head(gap[order(-gap$rho_gap_terms), ], 8), digits = 3, row.names = FALSE)

# Figure: rho per cell, by setting and term
s$cell_design <- ifelse(s$design %in% c("uniform_J100", "with_dist"), s$cell,
                        paste0(s$cell, " [", s$design, "]"))
s$setting <- factor(s$setting, c("oracle", "da_full", "da_seq"))
p <- ggplot(s, aes(x = rho, y = reorder(paste(block, cell_design), -as.integer(factor(block))),
                   colour = setting, shape = demonstration)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = rho - 1.96 * rho_mcse, xmax = rho + 1.96 * rho_mcse), width = 0, orientation = "y") +
  geom_point(size = 2) +
  facet_wrap(~ term) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c("test", "demonstration (expected to fail)")) +
  labs(x = "rho = estimate / target", y = NULL, colour = "Setting", shape = NULL,
       title = "Module A: recovery across group size, space and sampling",
       subtitle = "Target: own coefficient (oracle), pair surplus (deferred acceptance). 200 reps per cell, both ego sexes pooled.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "A_rho.pdf"), p, width = 11, height = 6.5)
cat("\nWrote A_summary.csv, A_checks.csv, A_term_gap.csv, A_rho.pdf\n")
