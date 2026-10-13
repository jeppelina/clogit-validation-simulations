# =============================================================================
# Module A (version 2): summary, pass criteria, figure. Run from the repo root
# after run_A.R. Writes results/v2_summary.csv, v2_checks.csv, v2_figure.pdf.
# rho = estimate / target; target = true taste (oracle) or pair surplus (2x,
# two-sided). For one_group, the majority's same-group term vs group k is
# -alt_gk, so its rho is -est / target.
# =============================================================================
library(ggplot2)
OUT <- "modules/A_recovery/results"
r <- do.call(rbind, lapply(file.path(OUT, paste0("v2_A", 0:3, ".csv")), read.csv))
b_true <- c(diag_g1 = 1, diag_g2 = 1, diag_g3 = 1, edu_diff = -1, ln_dist = -0.5,
            alt_g2 = -1, alt_g3 = -1)            # alt terms: one_group only, sign flipped below
r <- r[r$term %in% names(b_true), ]
r <- r[!(r$spec %in% c("full", "no_constants", "full_no_distance") & r$term %in% c("alt_g2", "alt_g3")), ]
mult <- ifelse(r$setting == "oracle", 1, 2)
r$rho <- r$est / (mult * b_true[r$term])
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
key <- c("block", "cell", "setting", "spec", "design", "term")
s <- do.call(rbind, lapply(split(r, r[, key], drop = TRUE), function(d)
  data.frame(d[1, key], rho = mean(d$rho, na.rm = TRUE), rho_median = median(d$rho, na.rm = TRUE),
             mcse = mcse(d$rho), n_ok = sum(!is.na(d$rho)), pair_rate = mean(d$pair_rate),
             mean_met = mean(d$mean_met), n_g3_egos = mean(d$n_g3_egos))))
rownames(s) <- NULL

# references: A0 full (same setting, term); loglinear (same cell, setting, term)
a0 <- s[s$block == "A0" & s$spec == "full", c("setting", "term", "rho", "mcse")]; names(a0)[3:4] <- c("rho_A0", "mcse_A0")
ll <- s[s$spec == "loglinear", c("cell", "setting", "term", "rho")]; names(ll)[4] <- "rho_loglinear"
s <- merge(s, a0, all.x = TRUE); s <- merge(s, ll, all.x = TRUE)
s$dev_A0 <- s$rho - s$rho_A0; s$dev_loglinear <- s$rho - s$rho_loglinear
s <- s[order(s$block, s$cell, match(s$setting, c("oracle", "gs", "seq")), s$spec, s$design, s$term), ]
write.csv(s, file.path(OUT, "v2_summary.csv"), row.names = FALSE)

# ---- pass criteria for `full` (with distance in A2; offset or uniform in A3) ----
tiny <- function(cell, term) cell == "90/8/2" & term == "diag_g3"
f <- s[s$spec == "full" & s$design != "strat_no_offset", ]
f$excluded <- tiny(f$cell, f$term)
chk <- list()
o <- f[f$setting == "oracle" & !f$excluded, ]
# Tolerance or 2.5 Monte Carlo SEs, whichever is larger (SPEC.md, amended
# before the production results were seen).
chk$oracle <- data.frame(check = "oracle |rho - 1| < max(0.02, 2.5 MCSE)", o[, c("block", "cell", "design", "term")],
                         value = o$rho - 1, pass = abs(o$rho - 1) < pmax(0.02, 2.5 * o$mcse))
t2 <- f[f$setting != "oracle" & grepl("diag", f$term) & !f$excluded, ]
eq <- do.call(rbind, lapply(split(t2, t2[, c("block", "cell", "setting", "design")], drop = TRUE), function(d)
  data.frame(check = "groups equal: range < max(0.05, 2.5 sqrt(2) max MCSE)", block = d$block[1], cell = d$cell[1],
             design = paste(d$setting[1], d$design[1]), term = "range of diag rho",
             value = diff(range(d$rho)), pass = diff(range(d$rho)) < max(0.05, 2.5 * sqrt(2) * max(d$mcse)))))
chk$equal <- eq
t3 <- f[f$setting != "oracle" & f$block != "A0" & !f$excluded & f$term != "ln_dist", ]
chk$stable <- data.frame(check = "within max(0.05, 2.5 MCSE) of A0", t3[, c("block", "cell", "design", "term")],
                         value = t3$dev_A0, pass = abs(t3$dev_A0) < pmax(0.05, 2.5 * sqrt(t3$mcse^2 + t3$mcse_A0^2)))
chk$stable$design <- paste(t3$setting, t3$design)
t4 <- f[f$setting != "oracle" & !is.na(f$dev_loglinear) & !f$excluded, ]
chk$ll <- data.frame(check = "agrees with log-linear within 0.05", t4[, c("block", "cell", "design", "term")],
                     value = t4$dev_loglinear, pass = abs(t4$dev_loglinear) < 0.05)
chk$ll$design <- paste(t4$setting, t4$design)
chk$oracle$design <- paste("oracle", chk$oracle$design)
checks <- do.call(rbind, lapply(chk, function(d) d[, c("check", "block", "cell", "design", "term", "value", "pass")]))
rownames(checks) <- NULL
write.csv(checks, file.path(OUT, "v2_checks.csv"), row.names = FALSE)

options(width = 200)
cat("\n== Pass rate by check ==\n")
print(aggregate(pass ~ check, checks, function(x) sprintf("%d / %d", sum(x), length(x))))
cat("\n== Failed checks ==\n")
print(checks[!checks$pass, ], digits = 3, row.names = FALSE)
cat("\n== full, same-group and education rho by cell ==\n")
w <- s[s$spec %in% c("full", "loglinear", "no_constants", "one_group"), c("block", "cell", "setting", "spec", "design", "term", "rho", "mcse", "n_ok")]
print(w, digits = 3, row.names = FALSE)

# ---- figure: same-group rho by group, full vs no_constants, over group shares ----
d <- s[s$block %in% c("A0", "A1") & s$spec %in% c("full", "no_constants") & grepl("diag", s$term) &
       !tiny(s$cell, s$term), ]
d$cell <- factor(d$cell, c("33/33/33", "70/20/10", "85/10/5", "90/8/2"))
d$group <- factor(d$term, paste0("diag_g", 1:3), c("largest group", "middle group", "smallest group"))
d$spec <- factor(d$spec, c("full", "no_constants"),
                 c("Well specified (alter-group constants)", "Without alter-group constants"))
d$setting <- factor(d$setting, c("oracle", "gs", "seq"), c("One-sided (oracle)", "Gale-Shapley", "Sequential search"))
p <- ggplot(d, aes(cell, rho, colour = group, group = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = rho - 1.96 * mcse, ymax = rho + 1.96 * mcse), width = 0) +
  facet_grid(spec ~ setting) +
  labs(x = "Group shares (largest / middle / smallest)", y = "rho = estimate / target", colour = NULL,
       title = "A well-specified conditional logit is not moved by group size",
       subtitle = "Identical tastes for every group. Target: true taste (one-sided), pair surplus (two-sided). 200 reps per cell.") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "v2_figure.pdf"), p, width = 10, height = 6)
cat("\nWrote v2_summary.csv, v2_checks.csv, v2_figure.pdf\n")
