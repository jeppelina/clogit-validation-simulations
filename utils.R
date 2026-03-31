# =============================================================================
# utils.R — Shared utility functions for choice model validation simulations
# =============================================================================
#
# This file provides the infrastructure used across all three modules:
#   - Agent generation (groups, traits, spatial placement)
#   - Workplace sorting (multi-stage DGP for Module 3)
#   - Choice set construction
#   - Partner choice simulation (based on random utility / Gumbel shocks)
#   - Model estimation wrapper around survival::clogit
#   - Parameter recovery assessment (bias, RMSE, coverage)
#   - Plotting helpers
#
# Usage: source("utils.R") at the top of each module script.
# =============================================================================

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

# =============================================================================
# 1. AGENT GENERATION
# =============================================================================

generate_agents <- function(N,
                            n_groups = 3,
                            group_proportions = c(0.7, 0.2, 0.1),
                            edu_means = c(12, 10, 14),
                            edu_sds = c(2, 3, 1.5),
                            include_openness = FALSE,
                            openness_edu_cor = 0.0) {
  #' Generate a population of agents with group membership and education.
  #'

  #' @param N Total population size
  #' @param n_groups Number of groups
  #' @param group_proportions Numeric vector of group proportions (must sum to 1)
  #' @param edu_means Mean education by group
  #' @param edu_sds SD of education by group
  #' @param include_openness If TRUE, add an unobserved "openness" trait

  #' @param openness_edu_cor Correlation between openness and education (Stage C)
  #' @return data.frame with columns: id, group, education, [openness]

  stopifnot(length(group_proportions) == n_groups)
  stopifnot(abs(sum(group_proportions) - 1) < 1e-6)

  # Assign groups
  group_sizes <- round(N * group_proportions)
  # Adjust for rounding
  group_sizes[1] <- N - sum(group_sizes[-1])

  agents <- data.frame(
    id = 1:N,
    group = rep(1:n_groups, times = group_sizes)
  )

  # Generate education by group
  agents$education <- NA_real_
  for (g in 1:n_groups) {
    idx <- agents$group == g
    n_g <- sum(idx)
    agents$education[idx] <- rnorm(n_g, mean = edu_means[g], sd = edu_sds[g])
  }

  # Optionally add unobserved openness trait
  if (include_openness) {
    if (abs(openness_edu_cor) > 0) {
      # Correlated with education
      # Use conditional normal: openness | education ~ N(rho * z_edu, sqrt(1 - rho^2))
      z_edu <- scale(agents$education)[, 1]
      agents$openness <- openness_edu_cor * z_edu +
        sqrt(1 - openness_edu_cor^2) * rnorm(N)
    } else {
      agents$openness <- rnorm(N)
    }
  }

  return(agents)
}


# =============================================================================
# 2. SPATIAL PLACEMENT (Module 2)
# =============================================================================

place_agents_on_grid <- function(agents, grid_size = 100,
                                 segregation_strength = 2.0) {
  #' Place agents on a 2D grid with trait-based residential clustering.
  #'
  #' Groups are assigned cluster centers. Agents are placed near their group's
  #' center(s) with noise controlled by segregation_strength.
  #'
  #' @param agents data.frame with id, group columns
  #' @param grid_size Size of the grid (grid_size x grid_size)
  #' @param segregation_strength Higher = more clustered (SD of placement noise)
  #' @return agents data.frame with added x, y columns

  n_groups <- length(unique(agents$group))

  # Assign group centers spread across the grid
  set.seed(NULL)  # ensure different centers each call if outer seed changes
  group_centers <- data.frame(
    group = 1:n_groups,
    cx = seq(grid_size * 0.2, grid_size * 0.8, length.out = n_groups),
    cy = rep(grid_size / 2, n_groups)
  )
  # Add some vertical spread for > 2 groups
  if (n_groups > 2) {
    group_centers$cy <- seq(grid_size * 0.3, grid_size * 0.7, length.out = n_groups)
  }

  # Place agents with noise around their group center
  spread <- grid_size / (segregation_strength * sqrt(n_groups))

  agents <- agents %>%
    left_join(group_centers, by = "group") %>%
    mutate(
      x = pmin(pmax(cx + rnorm(n(), 0, spread), 1), grid_size),
      y = pmin(pmax(cy + rnorm(n(), 0, spread), 1), grid_size)
    ) %>%
    select(-cx, -cy)

  return(agents)
}


# =============================================================================
# 3. WORKPLACE SORTING (Module 3)
# =============================================================================

sort_into_workplaces <- function(agents, n_workplaces = 20,
                                 edu_to_wp_strength = 2.0,
                                 openness_to_wp_strength = 0.0) {
  #' Sort agents into workplaces based on education and optionally openness.
  #'
  #' Workplaces have "education profiles" (centers on the education dimension).
  #' Agents are assigned probabilistically based on their distance to each
  #' workplace center. If openness_to_wp_strength > 0, workplaces also have
  #' openness profiles.
  #'
  #' @param agents data.frame with id, education, [openness]
  #' @param n_workplaces Number of workplaces
  #' @param edu_to_wp_strength How strongly education determines workplace
  #' @param openness_to_wp_strength How strongly openness determines workplace
  #' @return agents data.frame with added workplace column

  N <- nrow(agents)

  # Workplace centers on the education dimension
  edu_range <- range(agents$education)
  wp_centers_edu <- seq(edu_range[1] + 0.5, edu_range[2] - 0.5,
                        length.out = n_workplaces)

  # If using openness, create openness centers for workplaces
  if (openness_to_wp_strength > 0 && "openness" %in% names(agents)) {
    wp_centers_open <- rnorm(n_workplaces)
  }

  # Assign each agent to a workplace
  agents$workplace <- NA_integer_
  for (i in 1:N) {
    log_probs <- -edu_to_wp_strength * (agents$education[i] - wp_centers_edu)^2

    if (openness_to_wp_strength > 0 && "openness" %in% names(agents)) {
      log_probs <- log_probs -
        openness_to_wp_strength * (agents$openness[i] - wp_centers_open)^2
    }

    # Softmax
    log_probs <- log_probs - max(log_probs)  # numerical stability
    probs <- exp(log_probs) / sum(exp(log_probs))

    agents$workplace[i] <- sample(1:n_workplaces, 1, prob = probs)
  }

  return(agents)
}


# =============================================================================
# 4. CHOICE SET CONSTRUCTION
# =============================================================================

construct_choice_sets <- function(agents, J = 30, chooser_ids = NULL) {
  #' Construct choice sets by random sampling of alternatives.
  #'
  #' For each chooser, sample J alternatives from the population (excluding self).
  #' This mirrors the stratified random sampling approach in the empirical papers,
  #' simplified to uniform random sampling for the simulation.
  #'
  #' @param agents data.frame with id and all agent characteristics
  #' @param J Number of alternatives per choice set
  #' @param chooser_ids Optional: subset of agents who are choosers.
  #'   If NULL, all agents are choosers.
  #' @return data.frame in long format: one row per alternative per chooser,
  #'   with columns: chooser_id, alt_id, and all agent characteristics
  #'   prefixed with "ego_" (chooser) and "alt_" (alternative).

  if (is.null(chooser_ids)) {
    chooser_ids <- agents$id
  }

  N <- nrow(agents)
  all_ids <- agents$id

  # Pre-allocate list for speed
  cs_list <- vector("list", length(chooser_ids))

  for (idx in seq_along(chooser_ids)) {
    i <- chooser_ids[idx]

    # Sample alternatives (excluding self)
    available <- all_ids[all_ids != i]
    alts <- sample(available, size = min(J, length(available)), replace = FALSE)

    # Build choice set rows
    ego <- agents[agents$id == i, , drop = FALSE]
    alt_df <- agents[agents$id %in% alts, , drop = FALSE]

    cs <- data.frame(
      chooser_id = i,
      alt_id = alt_df$id
    )

    # Add ego characteristics
    for (col in setdiff(names(agents), "id")) {
      cs[[paste0("ego_", col)]] <- ego[[col]]
    }

    # Add alternative characteristics
    for (col in setdiff(names(agents), "id")) {
      cs[[paste0("alt_", col)]] <- alt_df[[col]]
    }

    cs_list[[idx]] <- cs
  }

  choice_data <- bind_rows(cs_list)
  return(choice_data)
}


# =============================================================================
# 5. COMPUTE MATCH VARIABLES
# =============================================================================

compute_match_variables <- function(choice_data) {
  #' Compute standard match-level variables from choice set data.
  #'
  #' Adds: edu_diff, same_group, same_workplace, [openness_diff], [distance]
  #' depending on what columns are available.

  cd <- choice_data

  # Education difference (absolute)
  if (all(c("ego_education", "alt_education") %in% names(cd))) {
    cd$edu_diff <- abs(cd$ego_education - cd$alt_education)
  }

  # Same group indicator
  if (all(c("ego_group", "alt_group") %in% names(cd))) {
    cd$same_group <- as.numeric(cd$ego_group == cd$alt_group)
  }

  # Same workplace indicator
  if (all(c("ego_workplace", "alt_workplace") %in% names(cd))) {
    cd$same_wp <- as.numeric(cd$ego_workplace == cd$alt_workplace)
  }

  # Openness difference (absolute)
  if (all(c("ego_openness", "alt_openness") %in% names(cd))) {
    cd$openness_diff <- abs(cd$ego_openness - cd$alt_openness)
  }

  # Euclidean distance (spatial)
  if (all(c("ego_x", "alt_x", "ego_y", "alt_y") %in% names(cd))) {
    cd$distance <- sqrt((cd$ego_x - cd$alt_x)^2 + (cd$ego_y - cd$alt_y)^2)
    cd$log_distance <- log(cd$distance + 1)
  }

  return(cd)
}


# =============================================================================
# 6. SIMULATE PARTNER CHOICE
# =============================================================================

simulate_choices <- function(choice_data, beta, var_names) {
  #' Simulate partner choices based on random utility model.
  #'
  #' For each choice set (chooser_id), compute utility for each alternative
  #' and select the one with highest utility (deterministic part + Gumbel shock).
  #'
  #' @param choice_data data.frame with chooser_id and match variables
  #' @param beta Named numeric vector of true coefficients
  #' @param var_names Character vector of variable names in choice_data
  #'   corresponding to the elements of beta
  #' @return choice_data with added "chosen" column (1/0)

  stopifnot(length(beta) == length(var_names))
  stopifnot(all(var_names %in% names(choice_data)))

  # Compute deterministic utility
  X <- as.matrix(choice_data[, var_names])
  V <- X %*% beta

  # Add Gumbel shocks
  n <- nrow(choice_data)
  eps <- -log(-log(runif(n)))  # standard Gumbel
  U <- V + eps

  # Select chosen alternative within each choice set
  choice_data$utility <- as.numeric(U)
  choice_data <- choice_data %>%
    group_by(chooser_id) %>%
    mutate(chosen = as.numeric(utility == max(utility))) %>%
    ungroup()

  # Handle ties (rare with continuous Gumbel, but just in case)
  # Keep only the first chosen if there are ties
  tie_groups <- choice_data %>%
    filter(chosen == 1) %>%
    group_by(chooser_id) %>%
    filter(n() > 1) %>%
    ungroup()

  if (nrow(tie_groups) > 0) {
    for (cid in unique(tie_groups$chooser_id)) {
      idx <- which(choice_data$chooser_id == cid & choice_data$chosen == 1)
      choice_data$chosen[idx[-1]] <- 0
    }
  }

  choice_data$utility <- NULL
  return(choice_data)
}


# =============================================================================
# 7. MODEL ESTIMATION WRAPPER
# =============================================================================

estimate_clogit <- function(choice_data, formula_str) {
  #' Estimate a conditional logit model using survival::clogit.
  #'
  #' @param choice_data data.frame with "chosen" and "chooser_id" columns
  #' @param formula_str Character string for the RHS of the formula
  #'   e.g., "edu_diff + same_group"
  #' @return List with: coefficients, se, z, p, converged, loglik, vcov, model

  formula <- as.formula(paste("chosen ~", formula_str, "+ strata(chooser_id)"))

  tryCatch({
    fit <- clogit(formula, data = choice_data, method = "efron")

    coefs <- coef(fit)
    se <- sqrt(diag(vcov(fit)))
    z <- coefs / se
    p <- 2 * pnorm(-abs(z))

    list(
      coefficients = coefs,
      se = se,
      z = z,
      p = p,
      converged = fit$info["convergence"] == 0,
      loglik = fit$loglik[2],
      vcov = vcov(fit),
      model = fit
    )
  }, error = function(e) {
    warning(paste("clogit failed:", e$message))
    list(
      coefficients = NA, se = NA, z = NA, p = NA,
      converged = FALSE, loglik = NA, vcov = NA, model = NULL
    )
  })
}


# =============================================================================
# 8. SINGLE SIMULATION RUN
# =============================================================================

run_single_sim <- function(sim_fn, ...) {
  #' Wrapper to run a single simulation and return results.
  #' sim_fn should accept ... and return a list with named coefficient vectors.
  sim_fn(...)
}


# =============================================================================
# 9. MONTE CARLO WRAPPER
# =============================================================================

run_monte_carlo <- function(sim_fn, R = 200, n_cores = 1, seed = 42, ...) {
  #' Run sim_fn R times and collect results.
  #'
  #' @param sim_fn Function that takes (...) and returns a list with
  #'   elements named by model (e.g., "M0", "M1"), each containing
  #'   a named numeric vector of estimated coefficients.
  #' @param R Number of Monte Carlo repetitions
  #' @param n_cores Number of cores for parallel execution
  #' @param seed Random seed (each rep gets seed + rep_id)
  #' @return List of data.frames, one per model, with columns:
  #'   rep, parameter, estimate

  if (n_cores > 1) {
    results <- mclapply(1:R, function(r) {
      set.seed(seed + r)
      tryCatch(sim_fn(...), error = function(e) NULL)
    }, mc.cores = n_cores)
  } else {
    results <- lapply(1:R, function(r) {
      set.seed(seed + r)
      tryCatch(sim_fn(...), error = function(e) NULL)
    })
  }

  # Remove failed runs
  results <- results[!sapply(results, is.null)]

  if (length(results) == 0) {
    stop("All simulation runs failed!")
  }

  cat(sprintf("Completed %d / %d runs successfully\n", length(results), R))

  return(results)
}


# =============================================================================
# 10. PARAMETER RECOVERY ASSESSMENT
# =============================================================================

assess_recovery <- function(mc_results, true_params, model_name = "M0") {
  #' Assess parameter recovery from Monte Carlo results.
  #'
  #' @param mc_results List of results from run_monte_carlo
  #' @param true_params Named numeric vector of true parameter values
  #' @param model_name Which model to assess (key in each result list)
  #' @return data.frame with: parameter, true, mean_est, bias, rel_bias_pct,
  #'   sd_est, rmse, coverage_95

  # Extract estimates
  param_names <- names(true_params)
  estimates <- matrix(NA, nrow = length(mc_results), ncol = length(param_names))
  ses <- matrix(NA, nrow = length(mc_results), ncol = length(param_names))

  for (r in seq_along(mc_results)) {
    res <- mc_results[[r]][[model_name]]
    if (!is.null(res) && !any(is.na(res$coefficients))) {
      for (p in seq_along(param_names)) {
        pname <- param_names[p]
        if (pname %in% names(res$coefficients)) {
          estimates[r, p] <- res$coefficients[pname]
          ses[r, p] <- res$se[pname]
        }
      }
    }
  }

  # Compute summary statistics
  recovery <- data.frame(
    parameter = param_names,
    true_value = as.numeric(true_params),
    mean_estimate = colMeans(estimates, na.rm = TRUE),
    median_estimate = apply(estimates, 2, median, na.rm = TRUE),
    sd_estimate = apply(estimates, 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  recovery$bias <- recovery$mean_estimate - recovery$true_value
  recovery$rel_bias_pct <- (recovery$bias / abs(recovery$true_value)) * 100
  recovery$rmse <- sqrt(recovery$bias^2 + recovery$sd_estimate^2)

  # Coverage: fraction of runs where true value is within 95% CI
  coverage <- numeric(length(param_names))
  for (p in seq_along(param_names)) {
    lower <- estimates[, p] - 1.96 * ses[, p]
    upper <- estimates[, p] + 1.96 * ses[, p]
    valid <- !is.na(lower)
    coverage[p] <- mean(true_params[p] >= lower[valid] &
                          true_params[p] <= upper[valid])
  }
  recovery$coverage_95 <- coverage

  recovery$n_valid <- colSums(!is.na(estimates))

  return(recovery)
}


# =============================================================================
# 11. MEDIATION DECOMPOSITION ASSESSMENT
# =============================================================================

assess_mediation <- function(mc_results, param_name = "edu_diff",
                             model_baseline = "M0", model_full = "M1") {
  #' Assess the mediation decomposition (coefficient change) across MC runs.
  #'
  #' @param mc_results List of results from run_monte_carlo
  #' @param param_name Name of the parameter whose change is the mediation
  #' @param model_baseline Baseline model name
  #' @param model_full Full model name
  #' @return data.frame with: mean_baseline, mean_full, mean_reduction,
  #'   mean_pct_reduction, sd_pct_reduction

  baseline_coefs <- numeric(length(mc_results))
  full_coefs <- numeric(length(mc_results))

  for (r in seq_along(mc_results)) {
    res <- mc_results[[r]]
    b0 <- res[[model_baseline]]$coefficients[param_name]
    b1 <- res[[model_full]]$coefficients[param_name]
    baseline_coefs[r] <- ifelse(is.null(b0) || is.na(b0), NA, b0)
    full_coefs[r] <- ifelse(is.null(b1) || is.na(b1), NA, b1)
  }

  valid <- !is.na(baseline_coefs) & !is.na(full_coefs)

  reduction <- baseline_coefs[valid] - full_coefs[valid]
  pct_reduction <- (reduction / baseline_coefs[valid]) * 100

  data.frame(
    mean_baseline_coef = mean(baseline_coefs[valid]),
    mean_full_coef = mean(full_coefs[valid]),
    mean_reduction = mean(reduction),
    mean_pct_reduction = mean(pct_reduction),
    median_pct_reduction = median(pct_reduction),
    sd_pct_reduction = sd(pct_reduction),
    n_valid = sum(valid)
  )
}


# =============================================================================
# 12. PLOTTING HELPERS
# =============================================================================

plot_parameter_recovery <- function(mc_results, true_params, model_name = "M0",
                                    title = "Parameter Recovery") {
  #' Plot distribution of estimated parameters vs true values.

  param_names <- names(true_params)
  estimates <- matrix(NA, nrow = length(mc_results), ncol = length(param_names))

  for (r in seq_along(mc_results)) {
    res <- mc_results[[r]][[model_name]]
    if (!is.null(res) && !any(is.na(res$coefficients))) {
      for (p in seq_along(param_names)) {
        pname <- param_names[p]
        if (pname %in% names(res$coefficients)) {
          estimates[r, p] <- res$coefficients[pname]
        }
      }
    }
  }

  df <- data.frame(estimates)
  names(df) <- param_names

  df_long <- df %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "estimate") %>%
    filter(!is.na(estimate))

  true_df <- data.frame(
    parameter = param_names,
    true_value = as.numeric(true_params)
  )

  ggplot(df_long, aes(x = estimate)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(data = true_df, aes(xintercept = true_value),
               color = "red", linewidth = 1, linetype = "dashed") +
    facet_wrap(~ parameter, scales = "free") +
    labs(title = title, x = "Estimated coefficient", y = "Count") +
    theme_minimal()
}


plot_bias_by_condition <- function(recovery_list, condition_var = "condition",
                                   title = "Bias by Condition") {
  #' Plot bias across experimental conditions.
  #'
  #' @param recovery_list List of data.frames from assess_recovery,
  #'   each with an added "condition" column.

  df <- bind_rows(recovery_list)

  ggplot(df, aes_string(x = condition_var, y = "bias", color = "parameter")) +
    geom_point(size = 3) +
    geom_line(aes(group = parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = bias - 1.96 * sd_estimate / sqrt(n_valid),
                      ymax = bias + 1.96 * sd_estimate / sqrt(n_valid)),
                  width = 0.1) +
    labs(title = title, x = condition_var, y = "Bias (estimate - true)") +
    theme_minimal()
}


# =============================================================================
# 13. DGP DIAGNOSTICS — Descriptive validation of each simulation step
# =============================================================================

describe_population <- function(agents) {
  #' Print descriptive statistics for the generated agent population.
  #' Shows education distributions overall and by group, plus openness if present.

  cat("\n  ---- POPULATION ----\n")
  cat(sprintf("  N = %d agents\n", nrow(agents)))

  # Education overall
  cat(sprintf("  Education: mean = %.2f, sd = %.2f, range = [%.1f, %.1f]\n",
              mean(agents$education), sd(agents$education),
              min(agents$education), max(agents$education)))

  # By group
  n_groups <- length(unique(agents$group))
  if (n_groups > 1) {
    cat(sprintf("  Groups: %d\n", n_groups))
    for (g in sort(unique(agents$group))) {
      sub <- agents[agents$group == g, ]
      cat(sprintf("    Group %d: n = %d (%.0f%%), edu mean = %.2f, sd = %.2f\n",
                  g, nrow(sub), 100 * nrow(sub) / nrow(agents),
                  mean(sub$education), sd(sub$education)))
    }
  } else {
    cat("  Groups: 1 (homogeneous population)\n")
  }

  # Openness
  if ("openness" %in% names(agents)) {
    cat(sprintf("  Openness: mean = %.2f, sd = %.2f\n",
                mean(agents$openness), sd(agents$openness)))
    rho <- cor(agents$education, agents$openness)
    cat(sprintf("  Cor(education, openness) = %.3f\n", rho))
  } else {
    cat("  Openness: not included\n")
  }
}


describe_workplaces <- function(agents) {
  #' Print descriptive statistics for workplace sorting.
  #' Shows size distribution, education segregation within workplaces.

  cat("\n  ---- WORKPLACES ----\n")

  wp_sizes <- table(agents$workplace)
  n_wp <- length(wp_sizes)
  cat(sprintf("  N workplaces: %d\n", n_wp))
  cat(sprintf("  Size distribution: min = %d, median = %d, max = %d\n",
              min(wp_sizes), median(wp_sizes), max(wp_sizes)))

  # Show the size of each workplace as a compact frequency row
  size_tbl <- sort(as.numeric(wp_sizes))
  cat(sprintf("  Sizes: %s\n", paste(size_tbl, collapse = ", ")))

  # Education within workplaces: mean within-workplace SD
  wp_edu_sd <- tapply(agents$education, agents$workplace, sd)
  wp_edu_sd <- wp_edu_sd[!is.na(wp_edu_sd)]  # drop singletons
  cat(sprintf("  Within-workplace edu SD: mean = %.2f (overall SD = %.2f)\n",
              mean(wp_edu_sd), sd(agents$education)))
  cat(sprintf("    -> Ratio = %.2f (lower = more segregated by education)\n",
              mean(wp_edu_sd) / sd(agents$education)))

  # Education means by workplace (sorted)
  wp_edu_mean <- sort(tapply(agents$education, agents$workplace, mean))
  cat(sprintf("  Workplace edu means range: [%.1f, %.1f]\n",
              min(wp_edu_mean), max(wp_edu_mean)))

  # If openness exists, show within-wp openness SD too
  if ("openness" %in% names(agents)) {
    wp_open_sd <- tapply(agents$openness, agents$workplace, sd)
    wp_open_sd <- wp_open_sd[!is.na(wp_open_sd)]
    cat(sprintf("  Within-workplace openness SD: mean = %.2f (overall SD = %.2f)\n",
                mean(wp_open_sd), sd(agents$openness)))

    # Correlation between wp-mean education and wp-mean openness
    wp_means <- data.frame(
      wp = names(tapply(agents$education, agents$workplace, mean)),
      edu = as.numeric(tapply(agents$education, agents$workplace, mean)),
      open = as.numeric(tapply(agents$openness, agents$workplace, mean))
    )
    cat(sprintf("  Cor(wp mean edu, wp mean openness) = %.3f\n",
                cor(wp_means$edu, wp_means$open)))
  }
}


describe_neighbourhoods <- function(agents) {
  #' Print descriptive statistics for neighbourhood assignment (if present).

  if (!"neighbourhood" %in% names(agents)) return(invisible(NULL))

  cat("\n  ---- NEIGHBOURHOODS ----\n")

  nbhd_sizes <- table(agents$neighbourhood)
  cat(sprintf("  N neighbourhoods: %d\n", length(nbhd_sizes)))
  cat(sprintf("  Size distribution: min = %d, median = %d, max = %d\n",
              min(nbhd_sizes), median(nbhd_sizes), max(nbhd_sizes)))

  nbhd_edu_sd <- tapply(agents$education, agents$neighbourhood, sd)
  nbhd_edu_sd <- nbhd_edu_sd[!is.na(nbhd_edu_sd)]
  cat(sprintf("  Within-neighbourhood edu SD: mean = %.2f (overall SD = %.2f)\n",
              mean(nbhd_edu_sd), sd(agents$education)))
}


describe_choice_sets <- function(choice_data) {
  #' Print descriptive statistics for the constructed choice sets.
  #' Shows choice set size, prevalence of key variables among alternatives.

  cat("\n  ---- CHOICE SETS ----\n")

  # Choice set sizes
  cs_sizes <- table(choice_data$chooser_id)
  cat(sprintf("  N choosers: %d, J per chooser: %d\n",
              length(cs_sizes), as.numeric(cs_sizes[1])))

  # Same workplace prevalence in choice sets
  if ("same_wp" %in% names(choice_data)) {
    wp_rate <- mean(choice_data$same_wp)
    cat(sprintf("  Pr(same workplace) in choice set: %.3f\n", wp_rate))
  }

  # Same neighbourhood prevalence
  if ("same_nbhd" %in% names(choice_data)) {
    nbhd_rate <- mean(choice_data$same_nbhd)
    cat(sprintf("  Pr(same neighbourhood) in choice set: %.3f\n", nbhd_rate))
  }

  # Same group prevalence
  if ("same_group" %in% names(choice_data)) {
    grp_rate <- mean(choice_data$same_group)
    cat(sprintf("  Pr(same group) in choice set: %.3f\n", grp_rate))
  }

  # Education difference distribution
  if ("edu_diff" %in% names(choice_data)) {
    cat(sprintf("  |edu_diff| in choice set: mean = %.2f, sd = %.2f\n",
                mean(choice_data$edu_diff), sd(choice_data$edu_diff)))
  }
}


describe_unions <- function(choice_data) {
  #' Print descriptive statistics for the simulated unions (chosen partners).
  #' Compares chosen alternatives to the full choice set baseline.

  cat("\n  ---- SIMULATED UNIONS ----\n")

  chosen <- choice_data[choice_data$chosen == 1, ]
  cat(sprintf("  N unions formed: %d\n", nrow(chosen)))

  # Same workplace: chosen vs baseline
  if ("same_wp" %in% names(choice_data)) {
    baseline_wp <- mean(choice_data$same_wp)
    chosen_wp <- mean(chosen$same_wp)
    cat(sprintf("  Same workplace:  chosen = %.3f  vs  baseline = %.3f  (ratio = %.1f)\n",
                chosen_wp, baseline_wp, chosen_wp / baseline_wp))
  }

  # Same neighbourhood: chosen vs baseline
  if ("same_nbhd" %in% names(choice_data)) {
    baseline_nbhd <- mean(choice_data$same_nbhd)
    chosen_nbhd <- mean(chosen$same_nbhd)
    cat(sprintf("  Same neighbourhood: chosen = %.3f  vs  baseline = %.3f  (ratio = %.1f)\n",
                chosen_nbhd, baseline_nbhd, chosen_nbhd / baseline_nbhd))
  }

  # Same group: chosen vs baseline
  if ("same_group" %in% names(choice_data)) {
    baseline_grp <- mean(choice_data$same_group)
    chosen_grp <- mean(chosen$same_group)
    cat(sprintf("  Same group:  chosen = %.3f  vs  baseline = %.3f  (ratio = %.1f)\n",
                chosen_grp, baseline_grp, chosen_grp / baseline_grp))
  }

  # Education difference: chosen vs baseline
  if ("edu_diff" %in% names(choice_data)) {
    baseline_edu <- mean(choice_data$edu_diff)
    chosen_edu <- mean(chosen$edu_diff)
    cat(sprintf("  |edu_diff|:  chosen = %.2f  vs  baseline = %.2f  (%.0f%% of baseline)\n",
                chosen_edu, baseline_edu, 100 * chosen_edu / baseline_edu))
  }

  # Openness difference if available
  if ("openness_diff" %in% names(choice_data)) {
    baseline_open <- mean(choice_data$openness_diff)
    chosen_open <- mean(chosen$openness_diff)
    cat(sprintf("  |openness_diff|:  chosen = %.2f  vs  baseline = %.2f\n",
                chosen_open, baseline_open))
  }
}


print_dgp_diagnostics <- function(agents, choice_data, label = "") {
  #' Run all descriptive checks for one DGP iteration.
  #' Call this once before the MC loop to show what the data looks like.

  cat(sprintf("\n  ====== DGP DIAGNOSTICS%s ======\n",
              if (nchar(label) > 0) paste0(": ", label) else ""))

  describe_population(agents)
  describe_workplaces(agents)
  describe_neighbourhoods(agents)
  describe_choice_sets(choice_data)
  describe_unions(choice_data)

  cat("  ====================================\n\n")
}


# =============================================================================
# 14. BILATERAL MATCHING — Two-sided matching with mutual acceptance
# =============================================================================
#
# In the one-sided DGP (Sections 4–6), each agent independently picks their
# best alternative from a random choice set. This is computationally clean
# but unrealistic: there is no mutual acceptance, no partner removal from
# the pool, and no competitive winnowing.
#
# The bilateral matching process adds realism:
#   1. All agents start unmatched
#   2. Each round: shuffle unmatched agents, pair them up sequentially
#   3. Both agents in a pair evaluate utility of the other
#   4. Both must have U(partner) > threshold to form a match
#   5. Matched pairs are removed from the pool
#   6. Repeat for T rounds
#   7. Remaining agents stay unmatched
#
# After matching, choice sets are constructed for clogit estimation:
#   - Each matched agent becomes a chooser
#   - Their actual partner is included in the choice set as "chosen"
#   - J-1 random alternatives are sampled from the full population
#
# IMPORTANT METHODOLOGICAL NOTES:
#   - Bilateral matching introduces selection: who matches depends on BOTH
#     partners' preferences. This means clogit with constructed choice sets
#     will show bias relative to the true DGP parameters.
#   - Python validation (30 MC reps) shows: beta_edu biased ~24% away from
#     zero (overestimated homophily), beta_same_group biased ~11% toward zero.
#   - This bias is INFORMATIVE: it mirrors the empirical setting where real
#     couples formed through bilateral selection, but we estimate one-sided
#     choice models with constructed choice sets.
#   - Status maximisation under bilateral random encounters does NOT produce
#     meaningful assortative matching (cor ≈ -0.04). Competitive winnowing
#     requires a full matching market (Gale-Shapley), not random encounters.
# =============================================================================


compute_pair_utility <- function(agents, ego_idx, alt_idx, beta, var_names) {
  #' Compute deterministic utility V of alt for ego.
  #'
  #' Supports the same variable names used in simulate_choices():
  #'   edu_diff, same_group, alt_edu, same_wp, openness_diff, distance,
  #'   log_distance, ego_x_alt_edu, same_nbhd
  #'
  #' @param agents data.frame with agent characteristics

  #' @param ego_idx Row index of the ego agent
  #' @param alt_idx Row index of the alternative agent
  #' @param beta Numeric vector of coefficients
  #' @param var_names Character vector of variable names (same length as beta)
  #' @return Scalar deterministic utility V

  V <- 0
  for (k in seq_along(var_names)) {
    val <- switch(var_names[k],
      "edu_diff" = abs(agents$education[ego_idx] - agents$education[alt_idx]),
      "same_group" = as.numeric(agents$group[ego_idx] == agents$group[alt_idx]),
      "alt_edu" = agents$education[alt_idx],
      "same_wp" = {
        if ("workplace" %in% names(agents))
          as.numeric(agents$workplace[ego_idx] == agents$workplace[alt_idx])
        else 0
      },
      "same_nbhd" = {
        if ("neighbourhood" %in% names(agents))
          as.numeric(agents$neighbourhood[ego_idx] == agents$neighbourhood[alt_idx])
        else 0
      },
      "openness_diff" = {
        if ("openness" %in% names(agents))
          abs(agents$openness[ego_idx] - agents$openness[alt_idx])
        else 0
      },
      "distance" = {
        if (all(c("x", "y") %in% names(agents)))
          sqrt((agents$x[ego_idx] - agents$x[alt_idx])^2 +
               (agents$y[ego_idx] - agents$y[alt_idx])^2)
        else 0
      },
      "log_distance" = {
        if (all(c("x", "y") %in% names(agents))) {
          d <- sqrt((agents$x[ego_idx] - agents$x[alt_idx])^2 +
                    (agents$y[ego_idx] - agents$y[alt_idx])^2)
          log(d + 1)
        } else 0
      },
      "ego_x_alt_edu" = agents$education[ego_idx] * agents$education[alt_idx],
      stop(paste("Unknown variable:", var_names[k]))
    )
    V <- V + beta[k] * val
  }
  return(V)
}


bilateral_matching <- function(agents, beta, var_names,
                               n_rounds = 20, threshold = -0.5,
                               encounter_fn = NULL,
                               diagnostics = FALSE,
                               verbose = FALSE) {
  #' Two-sided matching with bilateral acceptance.
  #'
  #' @param agents data.frame with id and agent characteristics
  #' @param beta Numeric vector of true preference coefficients
  #' @param var_names Character vector of match variable names
  #' @param n_rounds Number of encounter rounds
  #' @param threshold Minimum utility for acceptance (both must exceed)
  #'   Default -0.5 targets ~85-90% pairing with typical parameters.
  #' @param encounter_fn Optional function(unmatched_ids, agents) -> list of
  #'   pairs (2-column matrix). If NULL, uses random sequential pairing.
  #'   Use this for Module 2 to implement distance-based encounters.
  #' @param diagnostics If TRUE, store V_ab, V_ba (deterministic utilities)
  #'   for each matched pair. Used for pair-level analysis in Module 4.
  #' @param verbose Print progress during matching
  #' @return List with:
  #'   - matches: data.frame(ego_id, alt_id) of formed pairs
  #'   - unmatched: vector of unmatched agent ids
  #'   - pairing_rate: fraction of population that matched
  #'   - rounds_used: number of rounds actually run
  #'   - pair_diag: (if diagnostics=TRUE) data.frame with V_ab, V_ba per pair

  N <- nrow(agents)
  unmatched <- agents$id  # start with all agent IDs
  matches_ego <- integer(0)
  matches_alt <- integer(0)
  diag_V_ab <- numeric(0)
  diag_V_ba <- numeric(0)

  for (round_num in 1:n_rounds) {
    if (length(unmatched) < 2) break

    # --- Generate encounters ---
    if (is.null(encounter_fn)) {
      # Default: random sequential pairing
      shuffled <- sample(unmatched)
      n_pairs <- length(shuffled) %/% 2
      ego_ids <- shuffled[seq(1, 2 * n_pairs, by = 2)]
      alt_ids <- shuffled[seq(2, 2 * n_pairs, by = 2)]
      leftover <- if (length(shuffled) %% 2 == 1) shuffled[length(shuffled)] else integer(0)
    } else {
      # Custom encounter function (e.g., distance-based for Module 2)
      encounter_pairs <- encounter_fn(unmatched, agents)
      ego_ids <- encounter_pairs[, 1]
      alt_ids <- encounter_pairs[, 2]
      # Agents not in any pair stay unmatched
      paired_ids <- c(ego_ids, alt_ids)
      leftover <- setdiff(unmatched, paired_ids)
    }

    # --- Evaluate bilateral acceptance ---
    new_unmatched <- integer(0)

    for (p in seq_along(ego_ids)) {
      a <- ego_ids[p]
      b <- alt_ids[p]

      # Row indices in agents data.frame
      a_idx <- which(agents$id == a)
      b_idx <- which(agents$id == b)

      # Utility of b for a
      V_ab <- compute_pair_utility(agents, a_idx, b_idx, beta, var_names)
      eps_ab <- -log(-log(runif(1)))

      # Utility of a for b
      V_ba <- compute_pair_utility(agents, b_idx, a_idx, beta, var_names)
      eps_ba <- -log(-log(runif(1)))

      if (V_ab + eps_ab > threshold && V_ba + eps_ba > threshold) {
        matches_ego <- c(matches_ego, a)
        matches_alt <- c(matches_alt, b)
        if (diagnostics) {
          diag_V_ab <- c(diag_V_ab, V_ab)
          diag_V_ba <- c(diag_V_ba, V_ba)
        }
      } else {
        new_unmatched <- c(new_unmatched, a, b)
      }
    }

    unmatched <- c(new_unmatched, leftover)

    if (verbose && round_num %% 5 == 0) {
      cat(sprintf("    Round %d: %d pairs formed, %d still unmatched\n",
                  round_num, length(matches_ego), length(unmatched)))
    }
  }

  pairing_rate <- 2 * length(matches_ego) / N

  if (verbose) {
    cat(sprintf("  Bilateral matching complete: %d pairs (%.0f%%), %d unmatched\n",
                length(matches_ego), 100 * pairing_rate, length(unmatched)))
  }

  result <- list(
    matches = data.frame(ego_id = matches_ego, alt_id = matches_alt),
    unmatched = unmatched,
    pairing_rate = pairing_rate,
    rounds_used = min(n_rounds, which(length(unmatched) < 2)[1])
  )

  if (diagnostics && length(matches_ego) > 0) {
    result$pair_diag <- data.frame(
      ego_id = matches_ego,
      alt_id = matches_alt,
      V_ego_for_alt = diag_V_ab,  # deterministic utility of alt for ego
      V_alt_for_ego = diag_V_ba   # deterministic utility of ego for alt
    )
  }

  result
}


calibrate_threshold <- function(agents, beta, var_names,
                                target_acceptance_rate = 0.5,
                                n_sample_pairs = 5000, seed = 99) {
  #' Find the bilateral threshold that yields a target one-sided acceptance rate.
  #'

  #' For random pairs, computes U = V + Gumbel(0,1) and finds the threshold
  #' such that P(U > threshold) = target_acceptance_rate.
  #'
  #' The JOINT acceptance rate (both sides) is approximately target^2 for
  #' independent pairs, though correlated preferences make it slightly higher.
  #'
  #' @param agents data.frame with agent characteristics
  #' @param beta Numeric vector of true preference coefficients
  #' @param var_names Character vector of match variable names
  #' @param target_acceptance_rate Desired fraction of proposals accepted (one-sided)
  #' @param n_sample_pairs Number of random pairs to sample for calibration
  #' @param seed Random seed for reproducibility
  #' @return List with threshold, actual one-sided acceptance rate, V distribution summary

  set.seed(seed)
  N <- nrow(agents)
  n_sample_pairs <- min(n_sample_pairs, N * (N - 1) / 2)

  # Sample random pairs
  ego_idx <- sample(1:N, n_sample_pairs, replace = TRUE)
  alt_idx <- sample(1:N, n_sample_pairs, replace = TRUE)
  # Avoid self-pairs
  same <- ego_idx == alt_idx
  alt_idx[same] <- ((alt_idx[same]) %% N) + 1

  # Compute V for each pair
  V_vals <- numeric(n_sample_pairs)
  for (k in 1:n_sample_pairs) {
    V_vals[k] <- compute_pair_utility(agents, ego_idx[k], alt_idx[k], beta, var_names)
  }

  # Add Gumbel noise to get U = V + eps
  eps <- -log(-log(runif(n_sample_pairs)))
  U_vals <- V_vals + eps

  # Threshold = quantile of U such that P(U > threshold) = target
  # i.e., threshold = quantile(U, 1 - target)
  threshold <- as.numeric(quantile(U_vals, 1 - target_acceptance_rate))

  # Verify
  actual_rate <- mean(U_vals > threshold)

  list(
    threshold = threshold,
    target_acceptance_rate = target_acceptance_rate,
    actual_acceptance_rate = actual_rate,
    V_mean = mean(V_vals),
    V_sd = sd(V_vals),
    V_min = min(V_vals),
    V_max = max(V_vals)
  )
}


construct_choice_sets_bilateral <- function(agents, matches, J = 30) {
  #' Construct choice sets from bilateral matching results.
  #'
  #' For each matched pair (a, b):
  #'   - Agent a is a chooser with partner b as "chosen" + J-1 random alts
  #'   - Agent b is a chooser with partner a as "chosen" + J-1 random alts
  #'
  #' Unmatched agents are excluded (they have no observed choice).
  #'
  #' @param agents data.frame with id and all agent characteristics
  #' @param matches data.frame with ego_id and alt_id columns
  #' @param J Total choice set size (including the actual partner)
  #' @return data.frame in long format with chooser_id, alt_id, ego_/alt_
  #'   columns, and chosen indicator

  N <- nrow(agents)
  all_ids <- agents$id
  cs_list <- vector("list", 2 * nrow(matches))  # each pair -> 2 choice sets

  for (m in 1:nrow(matches)) {
    ego_a <- matches$ego_id[m]
    ego_b <- matches$alt_id[m]

    for (direction in 1:2) {
      # Direction 1: a chooses b; Direction 2: b chooses a
      if (direction == 1) {
        chooser <- ego_a; partner <- ego_b
      } else {
        chooser <- ego_b; partner <- ego_a
      }

      # Sample J-1 alternatives from population (excluding chooser and partner)
      available <- all_ids[all_ids != chooser & all_ids != partner]
      alts <- sample(available, size = min(J - 1, length(available)), replace = FALSE)
      alts <- c(alts, partner)  # partner goes last

      # Build choice set rows
      ego_row <- agents[agents$id == chooser, , drop = FALSE]
      alt_rows <- agents[match(alts, agents$id), , drop = FALSE]

      cs <- data.frame(chooser_id = chooser, alt_id = alts)

      # Add ego characteristics
      for (col in setdiff(names(agents), "id")) {
        cs[[paste0("ego_", col)]] <- ego_row[[col]]
      }

      # Add alternative characteristics
      for (col in setdiff(names(agents), "id")) {
        cs[[paste0("alt_", col)]] <- alt_rows[[col]]
      }

      # Mark the actual partner as chosen
      cs$chosen <- as.numeric(cs$alt_id == partner)

      list_idx <- (m - 1) * 2 + direction
      cs_list[[list_idx]] <- cs
    }
  }

  choice_data <- bind_rows(cs_list)
  return(choice_data)
}


construct_choice_sets_bilateral_from_pool <- function(agents, matches,
                                                       unmatched_pool = NULL,
                                                       J = 30) {
  #' Construct choice sets from bilateral matching, sampling alternatives

  #' from the unmatched pool rather than the full population.
  #'
  #' This is an alternative to construct_choice_sets_bilateral() that draws
  #' the J-1 non-chosen alternatives only from agents who remained unmatched.
  #' This more closely mirrors the "revealed preference" argument: the chooser
  #' could have matched with any of the unmatched agents but chose their
  #' actual partner instead.
  #'
  #' If unmatched_pool is NULL or has fewer than J-1 agents, falls back to
  #' sampling from the full population (same as construct_choice_sets_bilateral).
  #'
  #' @param agents data.frame with id and all agent characteristics
  #' @param matches data.frame with ego_id and alt_id columns
  #' @param unmatched_pool Vector of agent IDs who remained unmatched
  #' @param J Total choice set size (including the actual partner)
  #' @return data.frame in long format (same structure as construct_choice_sets_bilateral)

  N <- nrow(agents)
  all_ids <- agents$id
  cs_list <- vector("list", 2 * nrow(matches))

  # Determine the alternative pool
  use_unmatched <- !is.null(unmatched_pool) && length(unmatched_pool) >= (J - 1)

  for (m in 1:nrow(matches)) {
    ego_a <- matches$ego_id[m]
    ego_b <- matches$alt_id[m]

    for (direction in 1:2) {
      if (direction == 1) {
        chooser <- ego_a; partner <- ego_b
      } else {
        chooser <- ego_b; partner <- ego_a
      }

      # Sample J-1 alternatives from the specified pool
      if (use_unmatched) {
        pool <- unmatched_pool[unmatched_pool != chooser & unmatched_pool != partner]
      } else {
        pool <- all_ids[all_ids != chooser & all_ids != partner]
      }

      alts <- sample(pool, size = min(J - 1, length(pool)), replace = FALSE)
      alts <- c(alts, partner)  # partner goes last

      # Build choice set rows
      ego_row <- agents[agents$id == chooser, , drop = FALSE]
      alt_rows <- agents[match(alts, agents$id), , drop = FALSE]

      cs <- data.frame(chooser_id = chooser, alt_id = alts)

      for (col in setdiff(names(agents), "id")) {
        cs[[paste0("ego_", col)]] <- ego_row[[col]]
      }
      for (col in setdiff(names(agents), "id")) {
        cs[[paste0("alt_", col)]] <- alt_rows[[col]]
      }

      cs$chosen <- as.numeric(cs$alt_id == partner)

      list_idx <- (m - 1) * 2 + direction
      cs_list[[list_idx]] <- cs
    }
  }

  choice_data <- bind_rows(cs_list)
  return(choice_data)
}


# =============================================================================
# 15. DISTANCE-BASED ENCOUNTER FUNCTION (Module 2)
# =============================================================================

distance_encounter_fn <- function(max_encounter_distance = 30,
                                  distance_decay = 0.05) {
  #' Create a distance-based encounter function for spatial models.
  #'
  #' Returns a function suitable for the encounter_fn argument of
  #' bilateral_matching(). Agents are more likely to encounter
  #' nearby agents. The probability of encounter decays with distance.
  #'
  #' @param max_encounter_distance Maximum distance for an encounter
  #' @param distance_decay Decay rate for encounter probability
  #' @return A function(unmatched_ids, agents) -> 2-column matrix of pairs

  function(unmatched_ids, agents) {
    n <- length(unmatched_ids)
    if (n < 2) return(matrix(integer(0), ncol = 2))

    # Get positions of unmatched agents
    idx <- match(unmatched_ids, agents$id)
    x <- agents$x[idx]
    y <- agents$y[idx]

    # Shuffle and attempt to form distance-weighted pairs
    order <- sample(n)
    paired <- rep(FALSE, n)
    ego_ids <- integer(0)
    alt_ids <- integer(0)

    for (i in order) {
      if (paired[i]) next

      # Find unpaired agents within encounter distance
      candidates <- which(!paired & (1:n) != i)
      if (length(candidates) == 0) next

      dists <- sqrt((x[i] - x[candidates])^2 + (y[i] - y[candidates])^2)

      # Encounter probability decays with distance
      probs <- exp(-distance_decay * dists)
      probs[dists > max_encounter_distance] <- 0

      if (sum(probs) < 1e-10) next  # no one nearby

      probs <- probs / sum(probs)
      partner <- candidates[sample(length(candidates), 1, prob = probs)]

      ego_ids <- c(ego_ids, unmatched_ids[i])
      alt_ids <- c(alt_ids, unmatched_ids[partner])
      paired[i] <- TRUE
      paired[partner] <- TRUE
    }

    cbind(ego_ids, alt_ids)
  }
}


# =============================================================================
# 16. SEQUENTIAL ONE-SIDED MATCHING WITH PARTNER REMOVAL
# =============================================================================
#
# An intermediary DGP between one-sided (no removal) and bilateral (mutual
# acceptance). Agents choose sequentially: agent 1 picks their best
# alternative from J random candidates; the chosen partner is removed from
# the pool; agent 2 picks from the remaining pool; etc.
#
# This creates an opportunity structure effect: later choosers face a
# depleted pool (the "best" partners are already taken). But there is no
# bilateral acceptance — the chosen partner has no say.
# =============================================================================

sequential_onesided_matching <- function(agents, beta, var_names,
                                         J = 30, verbose = FALSE) {
  #' Sequential one-sided matching with partner removal.
  #'
  #' Agents are shuffled and each picks their best option from a random
  #' sample of J remaining agents. Chosen partners are removed.
  #'
  #' @param agents data.frame with id and agent characteristics
  #' @param beta Numeric vector of preference coefficients
  #' @param var_names Character vector of match variable names
  #' @param J Size of random choice set per chooser
  #' @param verbose Print progress
  #' @return Same format as bilateral_matching(): list with matches, unmatched, pairing_rate

  N <- nrow(agents)
  pool <- agents$id
  matches_ego <- integer(0)
  matches_alt <- integer(0)

  # Shuffle order of choosing
  choose_order <- sample(pool)

  for (chooser_id in choose_order) {
    # Remove self and already-matched agents from available pool
    available <- setdiff(pool, c(chooser_id, matches_ego, matches_alt))
    if (length(available) < 1) next

    # Sample J alternatives from available pool
    j_actual <- min(J, length(available))
    alt_ids <- sample(available, size = j_actual, replace = FALSE)

    # Compute utility of each alternative
    chooser_idx <- which(agents$id == chooser_id)
    utilities <- numeric(j_actual)
    for (k in seq_along(alt_ids)) {
      alt_idx <- which(agents$id == alt_ids[k])
      V <- compute_pair_utility(agents, chooser_idx, alt_idx, beta, var_names)
      eps <- -log(-log(runif(1)))
      utilities[k] <- V + eps
    }

    # Pick the best
    best <- which.max(utilities)
    matches_ego <- c(matches_ego, chooser_id)
    matches_alt <- c(matches_alt, alt_ids[best])
  }

  pairing_rate <- 2 * length(matches_ego) / N
  matched_ids <- c(matches_ego, matches_alt)
  unmatched <- setdiff(agents$id, matched_ids)

  if (verbose) {
    cat(sprintf("  Sequential one-sided: %d pairs (%.0f%%), %d unmatched\n",
                length(matches_ego), 100 * pairing_rate, length(unmatched)))
  }

  list(
    matches = data.frame(ego_id = matches_ego, alt_id = matches_alt),
    unmatched = unmatched,
    pairing_rate = pairing_rate
  )
}


# =============================================================================
# 17. DEFERRED ACCEPTANCE (GALE-SHAPLEY-LIKE) MATCHING
# =============================================================================
#
# The classical mechanism that produces stable matchings. Each round:
#   1. Each unmatched agent proposes to their most-preferred not-yet-rejected
#      partner from their preference list
#   2. Each agent who receives proposals tentatively accepts the best one
#      and rejects the rest
#   3. Tentatively matched agents can be displaced by better proposals
#
# This creates the competitive winnowing effect: high-value agents end up
# matched with each other because they keep displacing lower-value partners.
# Under status maximisation, this produces strong assortative matching.
#
# Simplified version: agents evaluate random samples (not full preference
# list) to make it computationally feasible with N=1000.
# =============================================================================

deferred_acceptance_matching <- function(agents, beta, var_names,
                                         J_proposals = 5,
                                         n_rounds = 30,
                                         verbose = FALSE) {
  #' Deferred acceptance (Gale-Shapley-like) matching.
  #'
  #' Each round, unmatched "proposers" propose to their best option from
  #' J_proposals random candidates. Targets tentatively accept the best
  #' proposal they've received (displacing any current tentative match).
  #' Displaced agents return to the proposer pool.
  #'
  #' @param agents data.frame with id and agent characteristics
  #' @param beta Numeric vector of preference coefficients
  #' @param var_names Character vector of match variable names
  #' @param J_proposals Number of candidates each proposer evaluates per round
  #' @param n_rounds Number of proposal rounds
  #' @param verbose Print progress
  #' @return Same format as bilateral_matching()

  N <- nrow(agents)

  # Track tentative matches: tentative_partner[i] = partner id or NA
  tentative_partner <- rep(NA_integer_, N)
  names(tentative_partner) <- agents$id

  # Track who has been proposed to (to avoid re-proposing)
  proposed_to <- vector("list", N)
  names(proposed_to) <- agents$id
  for (i in agents$id) proposed_to[[as.character(i)]] <- integer(0)

  for (round_num in 1:n_rounds) {
    # Find unmatched agents (proposers)
    unmatched_ids <- agents$id[is.na(tentative_partner)]
    if (length(unmatched_ids) < 2) break

    # Each unmatched agent proposes to their best available option
    proposals <- list()  # target_id -> list of proposer_ids

    for (proposer_id in unmatched_ids) {
      p_idx <- which(agents$id == proposer_id)
      already_proposed <- proposed_to[[as.character(proposer_id)]]

      # Available: everyone not yet proposed to, not self
      available <- setdiff(agents$id, c(proposer_id, already_proposed))
      if (length(available) < 1) next

      # Sample J_proposals candidates
      j_actual <- min(J_proposals, length(available))
      candidates <- sample(available, size = j_actual, replace = FALSE)

      # Evaluate and pick best
      utilities <- numeric(j_actual)
      for (k in seq_along(candidates)) {
        c_idx <- which(agents$id == candidates[k])
        V <- compute_pair_utility(agents, p_idx, c_idx, beta, var_names)
        eps <- -log(-log(runif(1)))
        utilities[k] <- V + eps
      }

      target <- candidates[which.max(utilities)]
      proposed_to[[as.character(proposer_id)]] <- c(already_proposed, target)

      # Record proposal
      target_key <- as.character(target)
      if (is.null(proposals[[target_key]])) {
        proposals[[target_key]] <- proposer_id
      } else {
        proposals[[target_key]] <- c(proposals[[target_key]], proposer_id)
      }
    }

    # Each target evaluates proposals (including current tentative partner)
    for (target_key in names(proposals)) {
      target_id <- as.integer(target_key)
      proposers <- proposals[[target_key]]
      t_idx <- which(agents$id == target_id)

      # Include current tentative partner in comparison
      current <- tentative_partner[target_key]
      all_candidates <- c(proposers, if (!is.na(current)) current)

      # Evaluate all candidates
      utilities <- numeric(length(all_candidates))
      for (k in seq_along(all_candidates)) {
        c_idx <- which(agents$id == all_candidates[k])
        V <- compute_pair_utility(agents, t_idx, c_idx, beta, var_names)
        eps <- -log(-log(runif(1)))
        utilities[k] <- V + eps
      }

      best <- all_candidates[which.max(utilities)]

      # Accept best, displace current if needed
      if (!is.na(current) && current != best) {
        # Displace current partner — they become unmatched again
        tentative_partner[as.character(current)] <- NA_integer_
      }

      # Reject all proposers except best (they stay unmatched)
      for (p in proposers) {
        if (p != best) {
          tentative_partner[as.character(p)] <- NA_integer_
        }
      }

      # Form tentative match
      tentative_partner[target_key] <- best
      tentative_partner[as.character(best)] <- target_id
    }

    if (verbose && round_num %% 10 == 0) {
      n_matched <- sum(!is.na(tentative_partner))
      cat(sprintf("    DA round %d: %d agents tentatively matched\n",
                  round_num, n_matched))
    }
  }

  # Convert tentative matches to pairs (avoid double-counting)
  matched_ids <- agents$id[!is.na(tentative_partner)]
  seen <- logical(N)
  names(seen) <- agents$id
  matches_ego <- integer(0)
  matches_alt <- integer(0)

  for (id in matched_ids) {
    if (seen[as.character(id)]) next
    partner <- tentative_partner[as.character(id)]
    if (!is.na(partner) && !seen[as.character(partner)]) {
      matches_ego <- c(matches_ego, id)
      matches_alt <- c(matches_alt, partner)
      seen[as.character(id)] <- TRUE
      seen[as.character(partner)] <- TRUE
    }
  }

  pairing_rate <- 2 * length(matches_ego) / N
  unmatched <- agents$id[is.na(tentative_partner)]

  if (verbose) {
    cat(sprintf("  Deferred acceptance: %d pairs (%.0f%%), %d unmatched\n",
                length(matches_ego), 100 * pairing_rate, length(unmatched)))
  }

  list(
    matches = data.frame(ego_id = matches_ego, alt_id = matches_alt),
    unmatched = unmatched,
    pairing_rate = pairing_rate
  )
}


# =============================================================================
# MESSAGE
# =============================================================================
cat("utils.R loaded successfully.\n")
cat("Available functions:\n")
cat("  generate_agents, place_agents_on_grid, sort_into_workplaces\n")
cat("  construct_choice_sets, compute_match_variables\n")
cat("  simulate_choices, estimate_clogit\n")
cat("  run_monte_carlo, assess_recovery, assess_mediation\n")
cat("  plot_parameter_recovery, plot_bias_by_condition\n")
cat("  describe_population, describe_workplaces, describe_choice_sets\n")
cat("  describe_unions, print_dgp_diagnostics\n")
cat("  compute_pair_utility, bilateral_matching\n")
cat("  construct_choice_sets_bilateral, construct_choice_sets_bilateral_from_pool\n")
cat("  distance_encounter_fn, sequential_onesided_matching, deferred_acceptance_matching\n")
