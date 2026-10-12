# =============================================================================
# Module 1: Group Size and Uneven Trait Distributions
# =============================================================================
#
# QUESTION: Can conditional logit recover correct homogamy and group-preference
# parameters when a continuous trait (education) is unevenly distributed across
# groups of different sizes?
#
# MOTIVATION: In the empirical data, education distributions differ sharply
# across ancestry groups. If the model conflates compositional effects with
# preference effects — especially for small groups with skewed trait
# distributions — the estimated ethnic preference coefficients could be biased.
#
# DGP:
#   - Multiple groups with different sizes (e.g., 70/15/15 or more extreme)
#   - Education distributed with different means and SDs across groups
#   - Partner choice utility:
#       U(i,j) = beta_edu * |edu_i - edu_j|
#              + beta_same_group * I(group_i == group_j)
#              + epsilon_ij  (Gumbel)
#
# EXPERIMENTS:
#   Exp 1.1: Balanced groups (33/33/33) — sanity check
#   Exp 1.2: Unbalanced groups (70/20/10) with same trait distribution
#   Exp 1.3: Unbalanced groups with different trait distributions
#   Exp 1.4: Extreme imbalance (85/10/5) with different trait distributions
#   Exp 1.5: Many small groups (50/10/10/10/5/5/5/5)
#
# For each experiment, we estimate:
#   Model A: edu_diff + same_group (the standard specification)
#   Model B: edu_diff + same_group + edu_diff:same_group (interaction)
#   Model C: edu_diff + group fixed effects (alt_group dummies)
#
# OUTPUT: Parameter recovery tables and bias plots across conditions.
# =============================================================================

source("utils.R")

# =============================================================================
# PARAMETERS
# =============================================================================

# Monte Carlo settings
R_SIMS <- 200        # number of repetitions per experiment
N_AGENTS <- 1000     # population size
J_ALTS <- 30         # alternatives per choice set
N_CORES <- 1         # set > 1 for parallel on multi-core machines

# True DGP parameters (same across experiments)
BETA_EDU <- -1.5       # homophily: prefer similar education
BETA_SAME_GROUP <- 0.8 # endogamy: prefer same group

# =============================================================================
# SINGLE SIMULATION FUNCTION
# =============================================================================

run_module1_sim <- function(group_proportions, edu_means, edu_sds,
                            N = N_AGENTS, J = J_ALTS,
                            beta_edu = BETA_EDU,
                            beta_same_group = BETA_SAME_GROUP,
                            matching_type = "one_sided",
                            bilateral_threshold = -0.5,
                            bilateral_rounds = 20) {
  #' Run one iteration of Module 1.
  #' Returns estimated coefficients for models A, B, C.
  #' @param matching_type "one_sided" (original) or "bilateral"

  n_groups <- length(group_proportions)

  # 1. Generate agents
  agents <- generate_agents(
    N = N,
    n_groups = n_groups,
    group_proportions = group_proportions,
    edu_means = edu_means,
    edu_sds = edu_sds
  )

  gen_beta <- c(beta_edu, beta_same_group)
  gen_vars <- c("edu_diff", "same_group")

  if (matching_type == "bilateral") {
    # Two-sided matching
    bm <- bilateral_matching(agents, gen_beta, gen_vars,
                             n_rounds = bilateral_rounds,
                             threshold = bilateral_threshold)
    if (nrow(bm$matches) < 20) return(NULL)
    choice_data <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    choice_data <- compute_match_variables(choice_data)
  } else {
    # One-sided (original)
    choice_data <- construct_choice_sets(agents, J = J)
    choice_data <- compute_match_variables(choice_data)
    choice_data <- simulate_choices(choice_data, beta = gen_beta, var_names = gen_vars)
  }

  # 5. Estimate models

  # Model A: Standard specification
  res_A <- estimate_clogit(choice_data, "edu_diff + same_group")

  # Model B: With interaction
  choice_data$edu_diff_x_same <- choice_data$edu_diff * choice_data$same_group
  res_B <- estimate_clogit(choice_data, "edu_diff + same_group + edu_diff_x_same")

  # Model C: Group fixed effects (alt_group dummies)
  # Create dummies for alternative's group (reference = group 1)
  for (g in 2:n_groups) {
    choice_data[[paste0("alt_group_", g)]] <- as.numeric(choice_data$alt_group == g)
  }
  group_fe_terms <- paste0("alt_group_", 2:n_groups, collapse = " + ")
  res_C <- estimate_clogit(choice_data,
                           paste("edu_diff +", group_fe_terms))

  list(
    model_A = res_A,
    model_B = res_B,
    model_C = res_C
  )
}


# =============================================================================
# EXPERIMENT DEFINITIONS
# =============================================================================

experiments <- list(
  # Exp 1.1: Balanced groups, same distributions
  exp1_1 = list(
    name = "1.1: Balanced, same distributions",
    group_proportions = c(1/3, 1/3, 1/3),
    edu_means = c(12, 12, 12),
    edu_sds = c(2, 2, 2)
  ),

  # Exp 1.2: Unbalanced groups, same distributions
  exp1_2 = list(
    name = "1.2: Unbalanced (70/20/10), same distributions",
    group_proportions = c(0.70, 0.20, 0.10),
    edu_means = c(12, 12, 12),
    edu_sds = c(2, 2, 2)
  ),

  # Exp 1.3: Unbalanced groups, different distributions
  exp1_3 = list(
    name = "1.3: Unbalanced (70/20/10), different distributions",
    group_proportions = c(0.70, 0.20, 0.10),
    edu_means = c(12, 10, 14),
    edu_sds = c(2, 3, 1.5)
  ),

  # Exp 1.4: Extreme imbalance, different distributions
  exp1_4 = list(
    name = "1.4: Extreme (85/10/5), different distributions",
    group_proportions = c(0.85, 0.10, 0.05),
    edu_means = c(12, 9, 15),
    edu_sds = c(2, 3, 1)
  ),

  # Exp 1.5: Many small groups
  exp1_5 = list(
    name = "1.5: Many groups (50/10/10/10/5/5/5/5)",
    group_proportions = c(0.50, 0.10, 0.10, 0.10, 0.05, 0.05, 0.05, 0.05),
    edu_means = c(12, 10, 14, 11, 9, 15, 13, 8),
    edu_sds = c(2, 3, 1.5, 2.5, 3, 1, 2, 3.5)
  )
)

# =============================================================================
# RUN EXPERIMENTS
# =============================================================================

cat("=" , rep("=", 59), "\n", sep = "")
cat("MODULE 1: Group Size and Uneven Trait Distributions\n")
cat("=" , rep("=", 59), "\n\n", sep = "")

all_results <- list()
all_recovery <- list()

for (exp_name in names(experiments)) {
  exp <- experiments[[exp_name]]
  cat(sprintf("\n--- Running %s ---\n", exp$name))
  cat(sprintf("    Group proportions: %s\n", paste(exp$group_proportions, collapse = ", ")))
  cat(sprintf("    Education means: %s\n", paste(exp$edu_means, collapse = ", ")))
  cat(sprintf("    Education SDs: %s\n", paste(exp$edu_sds, collapse = ", ")))

  # Run Monte Carlo
  mc_results <- run_monte_carlo(
    sim_fn = run_module1_sim,
    R = R_SIMS,
    n_cores = N_CORES,
    seed = 42,
    group_proportions = exp$group_proportions,
    edu_means = exp$edu_means,
    edu_sds = exp$edu_sds
  )

  all_results[[exp_name]] <- mc_results

  # Assess recovery for Model A (standard spec)
  true_A <- c(edu_diff = BETA_EDU, same_group = BETA_SAME_GROUP)
  rec_A <- assess_recovery(mc_results, true_A, model_name = "model_A")
  rec_A$experiment <- exp$name
  rec_A$model <- "Model A"

  cat("\nModel A (edu_diff + same_group):\n")
  print(rec_A[, c("parameter", "true_value", "mean_estimate", "bias",
                   "rel_bias_pct", "rmse", "coverage_95")])

  # Assess recovery for Model B (with interaction)
  # True interaction is 0 (DGP has no interaction)
  true_B <- c(edu_diff = BETA_EDU, same_group = BETA_SAME_GROUP, edu_diff_x_same = 0)
  rec_B <- assess_recovery(mc_results, true_B, model_name = "model_B")
  rec_B$experiment <- exp$name
  rec_B$model <- "Model B"

  cat("\nModel B (edu_diff + same_group + interaction):\n")
  print(rec_B[, c("parameter", "true_value", "mean_estimate", "bias",
                   "rel_bias_pct", "rmse", "coverage_95")])

  all_recovery[[paste0(exp_name, "_A")]] <- rec_A
  all_recovery[[paste0(exp_name, "_B")]] <- rec_B
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("SUMMARY: Bias in edu_diff across experiments (Model A)\n")
cat("=" , rep("=", 59), "\n\n", sep = "")

summary_df <- bind_rows(all_recovery) %>%
  filter(model == "Model A", parameter == "edu_diff") %>%
  select(experiment, mean_estimate, bias, rel_bias_pct, rmse, coverage_95)

print(summary_df)

cat("\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("SUMMARY: Bias in same_group across experiments (Model A)\n")
cat("=" , rep("=", 59), "\n\n", sep = "")

summary_df2 <- bind_rows(all_recovery) %>%
  filter(model == "Model A", parameter == "same_group") %>%
  select(experiment, mean_estimate, bias, rel_bias_pct, rmse, coverage_95)

print(summary_df2)

# =============================================================================
# PLOTS
# =============================================================================

# Save plots to results folder
dir.create("results", showWarnings = FALSE)

# Plot 1: Bias by experiment for Model A
recovery_df <- bind_rows(all_recovery) %>%
  filter(model == "Model A") %>%
  mutate(experiment = factor(experiment, levels = sapply(experiments, `[[`, "name")))

p1 <- ggplot(recovery_df, aes(x = experiment, y = bias, color = parameter)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = bias - 1.96 * sd_estimate / sqrt(n_valid),
                    ymax = bias + 1.96 * sd_estimate / sqrt(n_valid)),
                width = 0.2) +
  labs(title = "Module 1: Bias by Experimental Condition (Model A)",
       subtitle = paste0("True: beta_edu = ", BETA_EDU, ", beta_same_group = ", BETA_SAME_GROUP),
       x = "", y = "Bias (estimate - true)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

ggsave("results/module1_bias_by_condition.pdf", p1, width = 10, height = 6)
cat("\nPlot saved: results/module1_bias_by_condition.pdf\n")

# Plot 2: Coverage by experiment
p2 <- ggplot(recovery_df, aes(x = experiment, y = coverage_95, fill = parameter)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(title = "Module 1: 95% CI Coverage by Condition",
       x = "", y = "Coverage") +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

ggsave("results/module1_coverage.pdf", p2, width = 10, height = 6)
cat("Plot saved: results/module1_coverage.pdf\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

saveRDS(all_results, "results/module1_mc_results.rds")
saveRDS(bind_rows(all_recovery), "results/module1_recovery_summary.rds")
cat("\nResults saved to results/\n")

cat("\n--- Module 1 complete ---\n")
