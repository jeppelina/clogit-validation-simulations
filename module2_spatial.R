# =============================================================================
# Module 2: Spatial Opportunity vs. Homophily
# =============================================================================
#
# QUESTION: When education is spatially clustered and propinquity drives
# partner choice, does omitting propinquity from clogit create spurious
# educational homophily?
#
# HEADLINE: With global choice sets, spatial confounding is negligible.
# With local choice sets, it can be substantial.
#
# DGP:
#   - N agents on a 100×100 grid
#   - Education ~ N(12, 2), placed so that similar education → nearby location
#     (segregation_strength controls the tightness of spatial-education clustering)
#   - Partner choice utility:
#       U(i,j) = beta_homophily × |edu_i - edu_j|
#              + beta_propinquity × log(distance_ij + 1)
#              + epsilon_ij (Gumbel)
#   - Matching: oracle (one-sided)
#
# EXPERIMENTS:
#   2.1-2.3: Pure propinquity (no homophily), varying segregation
#            → does the reduced model show SPURIOUS homophily?
#   2.4a-c:  Both effects, varying segregation
#            → does omitting propinquity INFLATE homophily?
#   2.5:     Both effects + local choice sets (radius=20)
#            → does geographic restriction make it worse?
#   2.6:     Both effects + coarsened distance measurement
#            → does measurement error in distance leak into homophily?
#
# For each, estimate:
#   Full model:    edu_diff + log_distance
#   Reduced model: edu_diff  (omitting propinquity)
# =============================================================================

source("utils.R")

# =============================================================================
# PARAMETERS
# =============================================================================

R_SIMS <- 200
N_AGENTS <- 1000
J_ALTS <- 30
N_CORES <- 1
GRID_SIZE <- 100

# =============================================================================
# SPATIAL CHOICE SET CONSTRUCTION
# =============================================================================

construct_spatial_choice_sets <- function(agents, J = 30, radius = NULL) {
  #' Construct choice sets, optionally restricting to geographic radius.
  #'
  #' If radius is NULL, sample uniformly from all agents (standard).
  #' If radius is specified, only agents within that distance are eligible,
  #' and we sample from those (geographic dating market).

  N <- nrow(agents)
  cs_list <- vector("list", N)

  for (idx in 1:N) {
    i <- agents$id[idx]
    ego <- agents[agents$id == i, , drop = FALSE]

    if (!is.null(radius)) {
      # Distance to all other agents
      dists <- sqrt((agents$x - ego$x)^2 + (agents$y - ego$y)^2)
      available_mask <- (agents$id != i) & (dists <= radius)
      available <- agents$id[available_mask]

      if (length(available) < J) {
        # If not enough agents within radius, take all + sample rest
        outside <- agents$id[agents$id != i & !available_mask]
        n_needed <- J - length(available)
        extra <- sample(outside, size = min(n_needed, length(outside)), replace = FALSE)
        available <- c(available, extra)
      }
    } else {
      available <- agents$id[agents$id != i]
    }

    alts <- sample(available, size = min(J, length(available)), replace = FALSE)

    alt_df <- agents[agents$id %in% alts, , drop = FALSE]
    cs <- data.frame(chooser_id = i, alt_id = alt_df$id)

    for (col in setdiff(names(agents), "id")) {
      cs[[paste0("ego_", col)]] <- ego[[col]]
      cs[[paste0("alt_", col)]] <- alt_df[[col]]
    }

    cs_list[[idx]] <- cs
  }

  bind_rows(cs_list)
}


# =============================================================================
# SINGLE SIMULATION FUNCTIONS
# =============================================================================

run_module2_sim <- function(N = N_AGENTS, J = J_ALTS,
                            grid_size = GRID_SIZE,
                            segregation_strength = 2.0,
                            beta_homophily = -1.0,
                            beta_propinquity = -0.5,
                            beta_same_group = 0.0,
                            geographic_radius = NULL,
                            distance_noise_sd = 0,
                            matching_type = "one_sided",
                            bilateral_threshold = -0.5,
                            bilateral_rounds = 20,
                            bilateral_encounter_distance = 30,
                            bilateral_distance_decay = 0.05,
                            bilateral_random_encounters = FALSE) {
  #' Run one Module 2 iteration.
  #' @param matching_type "one_sided" or "bilateral"
  #' @param bilateral_encounter_distance Max distance for bilateral encounters
  #' @param bilateral_distance_decay Decay rate for encounter probability

  # 1. Generate agents (single group — spatial module focuses on education, not ethnicity)
  agents <- data.frame(
    id = 1:N,
    group = 1L,
    education = rnorm(N, mean = 12, sd = 2)
  )

  # 2. Place on grid with EDUCATION-BASED residential clustering
  #    Higher segregation_strength → stronger spatial sorting by education
  #    Location = education rank mapped to x-axis + noise controlled by segregation
  edu_rank <- rank(agents$education) / N  # 0-1 scale
  spread <- grid_size / (segregation_strength * 2)  # lower spread = tighter clustering
  agents$x <- pmin(pmax(edu_rank * grid_size + rnorm(N, 0, spread), 1), grid_size)
  agents$y <- runif(N, 1, grid_size)  # y is random (no structure)

  # Build DGP parameters
  gen_beta <- c(beta_homophily, beta_propinquity)
  gen_vars <- c("edu_diff", "log_distance")
  if (beta_same_group != 0) {
    gen_beta <- c(gen_beta, beta_same_group)
    gen_vars <- c(gen_vars, "same_group")
  }

  if (matching_type == "bilateral") {
    # --- Two-sided matching ---
    if (bilateral_random_encounters) {
      enc_fn <- NULL  # random encounters from full population
    } else {
      enc_fn <- distance_encounter_fn(
        max_encounter_distance = bilateral_encounter_distance,
        distance_decay = bilateral_distance_decay
      )
    }

    bm <- bilateral_matching(agents, gen_beta, gen_vars,
                             n_rounds = bilateral_rounds,
                             threshold = bilateral_threshold,
                             encounter_fn = enc_fn,
                             verbose = FALSE)
    if (nrow(bm$matches) < 20) return(NULL)

    cd <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    cd <- compute_match_variables(cd)

  } else {
    # --- One-sided (original) ---
    if (!is.null(geographic_radius)) {
      cd <- construct_spatial_choice_sets(agents, J = J, radius = geographic_radius)
    } else {
      cd <- construct_choice_sets(agents, J = J)
    }
    cd <- compute_match_variables(cd)

    # Add measurement error to distance if specified
    if (distance_noise_sd > 0) {
      cd$log_distance_true <- cd$log_distance
      cd$log_distance <- log(
        sqrt((round(cd$ego_x / distance_noise_sd) * distance_noise_sd -
              round(cd$alt_x / distance_noise_sd) * distance_noise_sd)^2 +
             (round(cd$ego_y / distance_noise_sd) * distance_noise_sd -
              round(cd$alt_y / distance_noise_sd) * distance_noise_sd)^2) + 1
      )
    }

    # Use true distance for DGP if we added noise
    if (distance_noise_sd > 0 && "log_distance_true" %in% names(cd)) {
      cd_gen <- cd
      cd_gen$log_distance <- cd_gen$log_distance_true
      cd_gen <- simulate_choices(cd_gen, beta = gen_beta, var_names = gen_vars)
      cd$chosen <- cd_gen$chosen
    } else {
      cd <- simulate_choices(cd, beta = gen_beta, var_names = gen_vars)
    }
  }

  # 6. Estimate models

  # 6a. Match statistics (chosen partners only)
  chosen <- cd[cd$chosen == 1, ]
  raw_dist <- sqrt((chosen$ego_x - chosen$alt_x)^2 + (chosen$ego_y - chosen$alt_y)^2)
  match_stats <- list(
    cor_ego_alt_edu = cor(chosen$ego_education, chosen$alt_education),
    mean_edu_diff   = mean(abs(chosen$ego_education - chosen$alt_education)),
    mean_distance   = mean(raw_dist),
    median_distance = median(raw_dist),
    spatial_edu_cor = cor(agents$education, agents$x)  # spatial clustering of education
  )

  # 6b. Full model: with propinquity
  if (beta_same_group != 0) {
    full_formula <- "edu_diff + log_distance + same_group"
  } else {
    full_formula <- "edu_diff + log_distance"
  }
  res_full <- estimate_clogit(cd, full_formula)

  # 6c. Reduced model: without propinquity
  if (beta_same_group != 0) {
    reduced_formula <- "edu_diff + same_group"
  } else {
    reduced_formula <- "edu_diff"
  }
  res_reduced <- estimate_clogit(cd, reduced_formula)

  list(full = res_full, reduced = res_reduced, match_stats = match_stats)
}


# =============================================================================
# EXPERIMENT DEFINITIONS
# =============================================================================

experiments <- list(
  # --- PURE PROPINQUITY (no homophily in DGP) ---
  # If segregation creates spatial-education correlation, and propinquity drives
  # choices, does the reduced model (omitting distance) show spurious homophily?

  exp2_1 = list(
    name = "2.1: Pure propinquity, low segregation",
    beta_homophily = 0, beta_propinquity = -0.8, beta_same_group = 0,
    segregation_strength = 1.0
  ),

  exp2_2 = list(
    name = "2.2: Pure propinquity, moderate segregation",
    beta_homophily = 0, beta_propinquity = -0.8, beta_same_group = 0,
    segregation_strength = 3.0
  ),

  exp2_3 = list(
    name = "2.3: Pure propinquity, high segregation",
    beta_homophily = 0, beta_propinquity = -0.8, beta_same_group = 0,
    segregation_strength = 6.0
  ),

  # --- BOTH EFFECTS (homophily + propinquity) ---
  # Does omitting propinquity inflate the edu_diff coefficient?

  exp2_4a = list(
    name = "2.4a: Both effects, low segregation",
    beta_homophily = -1.0, beta_propinquity = -0.5, beta_same_group = 0,
    segregation_strength = 1.0
  ),

  exp2_4b = list(
    name = "2.4b: Both effects, moderate segregation",
    beta_homophily = -1.0, beta_propinquity = -0.5, beta_same_group = 0,
    segregation_strength = 3.0
  ),

  exp2_4c = list(
    name = "2.4c: Both effects, high segregation",
    beta_homophily = -1.0, beta_propinquity = -0.5, beta_same_group = 0,
    segregation_strength = 6.0
  ),

  # --- GEOGRAPHIC CHOICE SET RESTRICTION ---
  # When choice sets are local (radius=20 on a 100x100 grid), does confounding bite?

  exp2_5 = list(
    name = "2.5: Both effects, high seg, local market (radius=20)",
    beta_homophily = -1.0, beta_propinquity = -0.5, beta_same_group = 0,
    segregation_strength = 6.0,
    geographic_radius = 20
  ),

  # --- DISTANCE MEASUREMENT ERROR ---
  exp2_6 = list(
    name = "2.6: Both effects, high seg, coarsened distance (grid=10)",
    beta_homophily = -1.0, beta_propinquity = -0.5, beta_same_group = 0,
    segregation_strength = 6.0,
    distance_noise_sd = 10
  )
)


# =============================================================================
# RUN EXPERIMENTS
# =============================================================================

cat("=" , rep("=", 59), "\n", sep = "")
cat("MODULE 2: Spatial Opportunity vs. Homophily\n")
cat("=" , rep("=", 59), "\n\n", sep = "")

dir.create("results", showWarnings = FALSE)

all_results <- list()
all_recovery <- list()

for (exp_name in names(experiments)) {
  exp <- experiments[[exp_name]]
  cat(sprintf("\n--- Running %s ---\n", exp$name))

  mc <- run_monte_carlo(
    sim_fn = run_module2_sim,
    R = R_SIMS,
    n_cores = N_CORES,
    seed = 42,
    segregation_strength = exp$segregation_strength,
    beta_homophily = exp$beta_homophily,
    beta_propinquity = exp$beta_propinquity,
    beta_same_group = exp$beta_same_group,
    geographic_radius = if (!is.null(exp$geographic_radius)) exp$geographic_radius else NULL,
    distance_noise_sd = if (!is.null(exp$distance_noise_sd)) exp$distance_noise_sd else 0
  )

  all_results[[exp_name]] <- mc

  # Assess full model
  true_full <- c(edu_diff = exp$beta_homophily, log_distance = exp$beta_propinquity)
  if (exp$beta_same_group != 0) {
    true_full <- c(true_full, same_group = exp$beta_same_group)
  }

  rec_full <- assess_recovery(mc, true_full, "full")
  rec_full$experiment <- exp$name
  rec_full$model <- "Full"

  # Assess reduced model (key test: does edu_diff absorb propinquity?)
  true_reduced <- c(edu_diff = exp$beta_homophily)
  if (exp$beta_same_group != 0) {
    true_reduced <- c(true_reduced, same_group = exp$beta_same_group)
  }
  rec_reduced <- assess_recovery(mc, true_reduced, "reduced")
  rec_reduced$experiment <- exp$name
  rec_reduced$model <- "Reduced"

  # Match statistics (across MC runs)
  valid_mc <- mc[sapply(mc, function(x) !is.null(x$match_stats))]
  if (length(valid_mc) > 0) {
    ms <- sapply(valid_mc, function(x) unlist(x$match_stats))
    cat("\n  Matching diagnostics (mean across MC runs):\n")
    cat(sprintf("    cor(ego_edu, alt_edu) = %.3f  (educational homogamy)\n",
                mean(ms["cor_ego_alt_edu", ])))
    cat(sprintf("    mean |edu_diff|       = %.2f\n", mean(ms["mean_edu_diff", ])))
    cat(sprintf("    mean partner distance = %.1f  (grid is %d×%d)\n",
                mean(ms["mean_distance", ]), GRID_SIZE, GRID_SIZE))
    cat(sprintf("    cor(education, x)     = %.3f  (spatial clustering of education)\n",
                mean(ms["spatial_edu_cor", ])))
  }

  cat("\n  Full model:\n")
  print(rec_full[, c("parameter", "true_value", "mean_estimate", "bias", "rel_bias_pct")])
  cat("\n  Reduced model (no propinquity):\n")
  print(rec_reduced[, c("parameter", "true_value", "mean_estimate", "bias", "rel_bias_pct")])

  all_recovery[[paste0(exp_name, "_full")]] <- rec_full
  all_recovery[[paste0(exp_name, "_red")]] <- rec_reduced

  gc(verbose = FALSE)
}


# =============================================================================
# BILATERAL MATCHING — "Both effects" experiments (2.4a-c) at 20% acceptance
# Tests whether spatial confounding interacts with bilateral selection bias.
# Encounters are random from the full population (not spatially structured).
# NOTE: Spatially structured encounters are a future extension.
# =============================================================================

cat("\n\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("BILATERAL MATCHING: Spatial × bilateral interaction\n")
cat("=" , rep("=", 59), "\n\n")

bilateral_exps <- c("exp2_4a", "exp2_4b", "exp2_4c")

for (exp_name in bilateral_exps) {
  exp <- experiments[[exp_name]]

  # Calibrate threshold for this DGP (edu_diff + log_distance)
  # Generate a calibration population with spatial placement
  set.seed(99)
  cal_agents <- data.frame(
    id = 1:N_AGENTS,
    group = 1L,
    education = rnorm(N_AGENTS, mean = 12, sd = 2)
  )
  edu_rank <- rank(cal_agents$education) / N_AGENTS
  spread <- GRID_SIZE / (exp$segregation_strength * 2)
  cal_agents$x <- pmin(pmax(edu_rank * GRID_SIZE + rnorm(N_AGENTS, 0, spread), 1), GRID_SIZE)
  cal_agents$y <- runif(N_AGENTS, 1, GRID_SIZE)

  cal_beta <- c(exp$beta_homophily, exp$beta_propinquity)
  cal_vars <- c("edu_diff", "log_distance")

  cal <- calibrate_threshold(cal_agents, cal_beta, cal_vars,
                             target_acceptance_rate = 0.20)
  bi_rounds <- min(ceiling(20 / 0.20^2), 500)
  rm(cal_agents); gc(verbose = FALSE)

  cat(sprintf("\n--- %s (bilateral 20%%, thr=%.3f) ---\n",
              exp$name, cal$threshold))

  mc_bi <- run_monte_carlo(
    sim_fn = run_module2_sim,
    R = R_SIMS,
    n_cores = N_CORES,
    seed = 42,
    segregation_strength = exp$segregation_strength,
    beta_homophily = exp$beta_homophily,
    beta_propinquity = exp$beta_propinquity,
    beta_same_group = exp$beta_same_group,
    matching_type = "bilateral",
    bilateral_threshold = cal$threshold,
    bilateral_rounds = bi_rounds,
    bilateral_random_encounters = TRUE
  )

  true_full <- c(edu_diff = exp$beta_homophily, log_distance = exp$beta_propinquity)

  rec_bi_full <- assess_recovery(mc_bi, true_full, "full")
  rec_bi_full$experiment <- paste(exp$name, "(bilateral 20%)")
  rec_bi_full$model <- "Full_bilateral"

  true_reduced <- c(edu_diff = exp$beta_homophily)
  rec_bi_red <- assess_recovery(mc_bi, true_reduced, "reduced")
  rec_bi_red$experiment <- paste(exp$name, "(bilateral 20%)")
  rec_bi_red$model <- "Reduced_bilateral"

  cat("\n  Full model (bilateral):\n")
  print(rec_bi_full[, c("parameter", "true_value", "mean_estimate", "bias", "rel_bias_pct")])
  cat("\n  Reduced model (bilateral):\n")
  print(rec_bi_red[, c("parameter", "true_value", "mean_estimate", "bias", "rel_bias_pct")])

  all_recovery[[paste0(exp_name, "_bi_full")]] <- rec_bi_full
  all_recovery[[paste0(exp_name, "_bi_red")]] <- rec_bi_red

  rm(mc_bi); gc(verbose = FALSE)
}


# =============================================================================
# KEY COMPARISON: edu_diff in full vs reduced across segregation levels
# =============================================================================

cat("\n\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("KEY RESULT: Homophily bias when propinquity is omitted\n")
cat("=" , rep("=", 59), "\n\n")

# Extract edu_diff estimates from experiments 2.4a-c (varying segregation, both effects)
seg_exps <- c("exp2_4a", "exp2_4b", "exp2_4c")
seg_labels <- c("Low", "Moderate", "High")

cat(sprintf("%-12s  %-15s  %-15s  %-15s\n",
            "Segregation", "Full (edu_diff)", "Reduced (edu_diff)", "Spurious bias"))
cat(rep("-", 65), "\n", sep = "")

for (i in seq_along(seg_exps)) {
  full_rec <- all_recovery[[paste0(seg_exps[i], "_full")]]
  red_rec <- all_recovery[[paste0(seg_exps[i], "_red")]]

  full_est <- full_rec$mean_estimate[full_rec$parameter == "edu_diff"]
  red_est <- red_rec$mean_estimate[red_rec$parameter == "edu_diff"]

  cat(sprintf("%-12s  %-15.3f  %-15.3f  %-15.3f\n",
              seg_labels[i], full_est, red_est, red_est - full_est))
}


# =============================================================================
# PLOTS
# =============================================================================

# Combine all recovery data
recovery_df <- bind_rows(all_recovery)

# Plot: Bias in edu_diff across experiments, full vs reduced
p_compare <- recovery_df %>%
  filter(parameter == "edu_diff") %>%
  ggplot(aes(x = experiment, y = bias, fill = model)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Module 2: Bias in Homophily Estimate",
       subtitle = "Full model (with propinquity) vs Reduced (without)",
       x = "", y = "Bias in edu_diff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave("results/module2_homophily_bias.pdf", p_compare, width = 10, height = 6)

# Plot: Propinquity recovery in full model
p_prop <- recovery_df %>%
  filter(parameter == "log_distance", model == "Full") %>%
  ggplot(aes(x = experiment, y = bias)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Module 2: Bias in Propinquity Estimate (Full Model)",
       x = "", y = "Bias in log_distance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave("results/module2_propinquity_bias.pdf", p_prop, width = 10, height = 6)

# =============================================================================
# SAVE
# =============================================================================

# Save CSV summary
write.csv(recovery_df, "results/module2_results.csv", row.names = FALSE)

cat("\nResults and plots saved to results/\n")
cat("\n--- Module 2 complete ---\n")
