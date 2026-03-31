# =============================================================================
# Module 3: Mediation via Sequential Opportunity Variables — Expanded
# =============================================================================
#
# QUESTION: Does the sequential addition of opportunity variables in clogit
# produce valid mediation estimates? Under what conditions does bias arise,
# and how large is it?
#
# This module tests the mediation decomposition by generating data under
# known DGPs and checking whether sequential model estimation recovers
# the true mediation share. It systematically varies both the confounding
# structure (Stages A–D) and a set of structural assumptions about the
# partner market.
#
# =============================================================================
# ASSUMPTIONS TESTED
# =============================================================================
#
# The current design makes many implicit choices. We identify each assumption,
# label it, and vary it systematically.
#
# --- ASSUMPTION 1: Mediation pathway strength ---
# How strongly education determines workplace sorting. Controls the true
# mediation share (how much of edu→partner runs through workplace).
# Values: 0.5 (weak), 2.0 (moderate, baseline), 5.0 (strong)
#
# --- ASSUMPTION 2: Choice set size (J) ---
# Number of alternatives sampled per chooser. Larger J provides more
# within-set variation but also changes the baseline probability and
# the effective "resolution" of the model. The empirical papers use
# varying J depending on market size.
# Values: 10 (small), 30 (baseline), 100 (large)
#
# --- ASSUMPTION 3: Confounder type ---
# How the unobserved trait enters the utility function.
# a) Homophily: beta * |openness_i - openness_j| (prefer similar)
# b) Level: beta * openness_j (prefer high-openness partners)
# c) Interaction: beta * openness_i * openness_j (assortative taste)
# Each type generates different patterns of confounding.
#
# --- ASSUMPTION 4: Group structure ---
# Whether agents come from multiple groups with different education
# distributions and workplace sorting rates. This mirrors the empirical
# setup where ethnic groups differ in education and institutional exposure.
# Values: single group (baseline) vs. 3 groups (70/20/10, different edu)
#
# --- ASSUMPTION 5: Multiple mediators ---
# Paper 4 adds variables sequentially: M0 → M1(+residence) → M2(+workplace)
# → M3(+university). With multiple mediators, confounding in one domain
# may propagate to the decomposition of others.
# Values: 1 mediator (baseline) vs. 2 mediators (workplace + neighbourhood)
#
# --- ASSUMPTION 6: Workplace concentration ---
# Whether workplaces are all similar size or highly unequal (a few large,
# many small). Affects the prevalence and variance of same_wp.
# Values: equal (20 wp, ~50 each) vs. concentrated (Zipf distribution)
#
# --- ASSUMPTION 7: Direct vs indirect effect ratio ---
# When the direct effect of education on partner choice (beta_edu) is
# large relative to the indirect pathway, the mediation share is small.
# The question: does confounding bias the mediation estimate more when
# the true share is small (inflating a small number) or when it's large?
# Varied by changing beta_edu while keeping the indirect pathway fixed.
#
# =============================================================================
# EXPERIMENTAL DESIGN
# =============================================================================
#
# We run two types of analyses:
#
# Part I:  CORE STAGES (A–D) with baseline assumptions
#          This replicates the original module with full MC.
#
# Part II: ASSUMPTION SWEEPS
#          For each assumption (1–7), we hold others at baseline and
#          vary the focal assumption across its levels. We run Stage B
#          (the core confounding case) for each level.
#
# Part III: INTERACTIONS
#          Selected cross-tabulations of assumptions where we expect
#          interactions (e.g., group structure × confounding type).
#
# =============================================================================
source("utils.R")

# =============================================================================
# GLOBAL PARAMETERS (BASELINE)
# =============================================================================

R_SIMS <- 200       # MC reps per condition
N_AGENTS <- 1000
J_ALTS <- 30
N_CORES <- 1
N_WORKPLACES <- 20

# True DGP parameters (baseline)
BETA_EDU <- -1.0          # direct effect of education similarity
BETA_WP <- 1.5            # effect of shared workplace
BETA_OPENNESS <- -0.8     # unobserved confounder effect (homophily type)
EDU_TO_WP <- 2.0          # education → workplace sorting strength
OPENNESS_TO_WP <- 1.0     # confounding: openness → workplace
OPENNESS_EDU_COR <- 0.0   # correlation between openness and education


# =============================================================================
# HELPER: Generate workplaces with unequal sizes (Assumption 6)
# =============================================================================

sort_into_workplaces_concentrated <- function(agents, n_workplaces = 20,
                                               edu_to_wp_strength = 2.0,
                                               openness_to_wp_strength = 0.0,
                                               concentration = "equal") {
  #' Sort agents into workplaces with optionally concentrated sizes.
  #'
  #' concentration = "equal": all workplaces equally attractive (baseline)
  #' concentration = "zipf": workplaces have Zipf-distributed attractiveness
  #'   (a few large workplaces, many small ones)

  N <- nrow(agents)
  edu_range <- range(agents$education)
  wp_centers_edu <- seq(edu_range[1] + 0.5, edu_range[2] - 0.5,
                        length.out = n_workplaces)

  # Workplace size weights
  if (concentration == "zipf") {
    # Zipf: rank^(-1) distribution
    wp_size_weight <- log(1 / (1:n_workplaces))  # log scale for softmax
  } else {
    wp_size_weight <- rep(0, n_workplaces)
  }

  if (openness_to_wp_strength > 0 && "openness" %in% names(agents)) {
    wp_centers_open <- rnorm(n_workplaces)
  }

  agents$workplace <- NA_integer_
  for (i in 1:N) {
    log_probs <- -edu_to_wp_strength * (agents$education[i] - wp_centers_edu)^2 +
      wp_size_weight

    if (openness_to_wp_strength > 0 && "openness" %in% names(agents)) {
      log_probs <- log_probs -
        openness_to_wp_strength * (agents$openness[i] - wp_centers_open)^2
    }

    log_probs <- log_probs - max(log_probs)
    probs <- exp(log_probs) / sum(exp(log_probs))
    agents$workplace[i] <- sample(1:n_workplaces, 1, prob = probs)
  }

  return(agents)
}


# =============================================================================
# HELPER: Assign neighbourhoods (Assumption 5 — multiple mediators)
# =============================================================================

assign_neighbourhoods <- function(agents, n_nbhds = 10,
                                  edu_to_nbhd = 1.5,
                                  openness_to_nbhd = 0.0) {
  #' Assign agents to neighbourhoods based on education (and optionally openness).
  #' Neighbourhoods are conceptually residential areas with education profiles.

  N <- nrow(agents)
  edu_range <- range(agents$education)
  nbhd_centers <- seq(edu_range[1] + 0.5, edu_range[2] - 0.5,
                       length.out = n_nbhds)

  if (openness_to_nbhd > 0 && "openness" %in% names(agents)) {
    nbhd_open_centers <- rnorm(n_nbhds)
  }

  agents$neighbourhood <- NA_integer_
  for (i in 1:N) {
    log_probs <- -edu_to_nbhd * (agents$education[i] - nbhd_centers)^2

    if (openness_to_nbhd > 0 && "openness" %in% names(agents)) {
      log_probs <- log_probs -
        openness_to_nbhd * (agents$openness[i] - nbhd_open_centers)^2
    }

    log_probs <- log_probs - max(log_probs)
    probs <- exp(log_probs) / sum(exp(log_probs))
    agents$neighbourhood[i] <- sample(1:n_nbhds, 1, prob = probs)
  }

  return(agents)
}


# =============================================================================
# CORE SIMULATION FUNCTION (handles all assumption variations)
# =============================================================================

run_module3_sim <- function(
    # Population
    N = N_AGENTS,
    J = J_ALTS,
    n_workplaces = N_WORKPLACES,

    # True parameters
    beta_edu = BETA_EDU,
    beta_wp = BETA_WP,
    beta_openness = BETA_OPENNESS,

    # Mediation pathway (A1)
    edu_to_wp = EDU_TO_WP,

    # Confounding
    openness_to_wp = OPENNESS_TO_WP,
    openness_edu_cor = OPENNESS_EDU_COR,
    include_confounder = TRUE,

    # Confounder type (A3)
    confounder_type = "homophily",  # "homophily", "level", "interaction"

    # Group structure (A4)
    use_groups = FALSE,
    group_proportions = c(0.7, 0.2, 0.1),
    edu_means = c(12, 10, 14),
    edu_sds = c(2, 3, 1.5),

    # Multiple mediators (A5)
    use_neighbourhood = FALSE,
    beta_nbhd = 1.0,
    edu_to_nbhd = 1.5,
    openness_to_nbhd = 0.0,

    # Workplace concentration (A6)
    wp_concentration = "equal",

    # Choice set selection (Stage D)
    openness_cs_strength = 0.0,

    # Matching type (one-sided vs bilateral)
    matching_type = "one_sided",       # "one_sided" or "bilateral"
    bilateral_threshold = -0.5,
    bilateral_rounds = 20,

    # Diagnostics
    verbose = FALSE
) {
  #' Run one Module 3 iteration with flexible assumptions.
  #' @param matching_type "one_sided" (original) or "bilateral"
  #' @param bilateral_threshold Minimum utility for bilateral acceptance
  #' @param bilateral_rounds Number of bilateral encounter rounds

  # --- 1. GENERATE AGENTS ---
  if (use_groups) {
    n_groups <- length(group_proportions)
    agents <- generate_agents(
      N = N, n_groups = n_groups,
      group_proportions = group_proportions,
      edu_means = edu_means, edu_sds = edu_sds,
      include_openness = include_confounder,
      openness_edu_cor = openness_edu_cor
    )
  } else {
    agents <- generate_agents(
      N = N, n_groups = 1,
      group_proportions = 1,
      edu_means = 12, edu_sds = 2,
      include_openness = include_confounder,
      openness_edu_cor = openness_edu_cor
    )
  }

  # --- 2. SORT INTO WORKPLACES ---
  agents <- sort_into_workplaces_concentrated(
    agents,
    n_workplaces = n_workplaces,
    edu_to_wp_strength = edu_to_wp,
    openness_to_wp_strength = if (include_confounder) openness_to_wp else 0,
    concentration = wp_concentration
  )

  # --- 3. OPTIONALLY ASSIGN NEIGHBOURHOODS ---
  if (use_neighbourhood) {
    agents <- assign_neighbourhoods(
      agents, n_nbhds = 10,
      edu_to_nbhd = edu_to_nbhd,
      openness_to_nbhd = if (include_confounder) openness_to_nbhd else 0
    )
  }

  # --- 4. BUILD DGP PARAMETER VECTOR ---
  # (needed for both one-sided simulation and bilateral matching)
  gen_beta <- c(beta_edu)
  gen_vars <- c("edu_diff")

  gen_beta <- c(gen_beta, beta_wp)
  gen_vars <- c(gen_vars, "same_wp")

  if (use_neighbourhood) {
    gen_beta <- c(gen_beta, beta_nbhd)
    gen_vars <- c(gen_vars, "same_nbhd")
  }

  # Determine confounder variable type and coefficient
  if (include_confounder && "openness" %in% names(agents)) {
    if (confounder_type == "homophily") {
      open_beta <- beta_openness  # negative = prefer similar
      open_var <- "openness_diff"
    } else if (confounder_type == "level") {
      open_beta <- abs(beta_openness)  # positive = prefer high-openness
      open_var <- "alt_openness"  # handled specially below for bilateral
    } else if (confounder_type == "interaction") {
      open_beta <- abs(beta_openness)  # positive = assortative matching
      open_var <- "ego_x_alt_openness"  # handled specially below
    }
  }


  # --- 5. MATCHING TYPE BRANCH ---
  if (matching_type == "bilateral") {
    # =====================================================================
    # BILATERAL MATCHING PATH
    # =====================================================================
    # For bilateral matching, the confounder enters compute_pair_utility().
    # We need to use variable names that compute_pair_utility() understands.
    # "openness_diff" is already supported. For "level" and "interaction"
    # confounder types, we need to handle them in the choice set construction.

    bilateral_gen_beta <- gen_beta
    bilateral_gen_vars <- gen_vars

    if (include_confounder && "openness" %in% names(agents)) {
      if (confounder_type == "homophily") {
        bilateral_gen_beta <- c(bilateral_gen_beta, open_beta)
        bilateral_gen_vars <- c(bilateral_gen_vars, "openness_diff")
      }
      # For "level" and "interaction" types, we skip them in bilateral matching
      # utility (they operate differently and are handled post-matching).
      # This is a simplification — the bilateral matching uses the primary
      # preference variables (edu_diff, same_wp, same_nbhd) + openness homophily
      # only. The choice-set-based estimation still tests all confounder types.
    }

    bm <- bilateral_matching(agents, bilateral_gen_beta, bilateral_gen_vars,
                             n_rounds = bilateral_rounds,
                             threshold = bilateral_threshold,
                             verbose = verbose)

    if (nrow(bm$matches) < 20) {
      # Not enough matches — return NA results
      return(list(M0 = list(coefficients = NA, se = NA),
                  M1 = list(coefficients = NA, se = NA)))
    }

    # Construct choice sets from bilateral matching outcomes
    cd <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    cd <- compute_match_variables(cd)

    # Add neighbourhood variable if applicable
    if (use_neighbourhood && "ego_neighbourhood" %in% names(cd)) {
      cd$same_nbhd <- as.numeric(cd$ego_neighbourhood == cd$alt_neighbourhood)
    }

    # Add confounder variable for estimation (not for matching — already done)
    if (include_confounder && "ego_openness" %in% names(cd)) {
      if (confounder_type == "homophily") {
        cd$openness_var <- abs(cd$ego_openness - cd$alt_openness)
      } else if (confounder_type == "level") {
        cd$openness_var <- cd$alt_openness
      } else if (confounder_type == "interaction") {
        cd$openness_var <- cd$ego_openness * cd$alt_openness
      }
    }

    if (verbose) {
      cat(sprintf("  Bilateral matching: %d pairs, rate = %.0f%%\n",
                  nrow(bm$matches), 100 * bm$pairing_rate))
    }

  } else {
    # =====================================================================
    # ONE-SIDED PATH (original)
    # =====================================================================
    if (openness_cs_strength > 0 && "openness" %in% names(agents)) {
      # Custom construction with openness-based selection
      all_ids <- agents$id
      log_weights <- openness_cs_strength * agents$openness
      log_weights <- log_weights - max(log_weights)
      base_weights <- exp(log_weights)

      cs_list <- vector("list", N)
      for (idx in 1:N) {
        i <- agents$id[idx]
        available <- all_ids[all_ids != i]
        w <- base_weights[all_ids != i]
        w <- w / sum(w)

        alts <- sample(available, size = min(J, length(available)),
                       replace = FALSE, prob = w)

        ego <- agents[agents$id == i, , drop = FALSE]
        alt_df <- agents[agents$id %in% alts, , drop = FALSE]

        cs <- data.frame(chooser_id = i, alt_id = alt_df$id)
        for (col in setdiff(names(agents), "id")) {
          cs[[paste0("ego_", col)]] <- ego[[col]]
          cs[[paste0("alt_", col)]] <- alt_df[[col]]
        }
        cs_list[[idx]] <- cs
      }
      cd <- bind_rows(cs_list)
    } else {
      cd <- construct_choice_sets(agents, J = J)
    }

    cd <- compute_match_variables(cd)

    # Add neighbourhood variable if applicable
    if (use_neighbourhood) {
      cd$same_nbhd <- as.numeric(cd$ego_neighbourhood == cd$alt_neighbourhood)
    }

    # Compute confounder variable (for DGP, not observed)
    if (include_confounder && "ego_openness" %in% names(cd)) {
      if (confounder_type == "homophily") {
        cd$openness_var <- abs(cd$ego_openness - cd$alt_openness)
      } else if (confounder_type == "level") {
        cd$openness_var <- cd$alt_openness
      } else if (confounder_type == "interaction") {
        cd$openness_var <- cd$ego_openness * cd$alt_openness
      }
    }

    # Build simulation beta vector including confounder
    sim_beta <- gen_beta
    sim_vars <- gen_vars
    if (include_confounder && "openness_var" %in% names(cd)) {
      sim_beta <- c(sim_beta, open_beta)
      sim_vars <- c(sim_vars, "openness_var")
    }

    cd <- simulate_choices(cd, beta = sim_beta, var_names = sim_vars)

    if (verbose) {
      print_dgp_diagnostics(agents, cd)
    }
  }


  # --- 6. ESTIMATE SEQUENTIAL MODELS ---
  # M0: education only
  M0 <- estimate_clogit(cd, "edu_diff")

  # M1: education + workplace
  M1 <- estimate_clogit(cd, "edu_diff + same_wp")

  results <- list(M0 = M0, M1 = M1)

  # If multiple mediators: M2 = edu + workplace + neighbourhood
  if (use_neighbourhood && "same_nbhd" %in% names(cd)) {
    M1b <- estimate_clogit(cd, "edu_diff + same_nbhd")
    M2 <- estimate_clogit(cd, "edu_diff + same_wp + same_nbhd")
    results$M1b <- M1b  # neighbourhood only
    results$M2 <- M2    # both mediators
  }

  # If groups: also estimate with group controls
  if (use_groups) {
    n_groups <- length(unique(agents$group))
    for (g in 2:n_groups) {
      cd[[paste0("alt_group_", g)]] <- as.numeric(cd$alt_group == g)
    }
    group_terms <- paste0("alt_group_", 2:n_groups, collapse = " + ")

    M0g <- estimate_clogit(cd, paste("edu_diff +", group_terms))
    M1g <- estimate_clogit(cd, paste("edu_diff + same_wp +", group_terms))
    results$M0g <- M0g
    results$M1g <- M1g
  }

  return(results)
}


# =============================================================================
# PART I: CORE STAGES (A–D) WITH BASELINE ASSUMPTIONS
# =============================================================================

run_part_I <- function() {
  cat("\n")
  cat("=" , rep("=", 70), "\n", sep = "")
  cat("PART I: Core Stages A-D (baseline assumptions)\n")
  cat("=" , rep("=", 70), "\n")

  true_M1 <- c(edu_diff = BETA_EDU, same_wp = BETA_WP)

  all_results <- list()

  # Calibrate bilateral threshold at 20% acceptance rate for this DGP
  set.seed(99)
  cal_agents <- generate_agents(N = N_AGENTS, n_groups = 1,
                                group_proportions = 1,
                                edu_means = 12, edu_sds = 2,
                                include_openness = TRUE)
  cal_agents <- sort_into_workplaces_concentrated(cal_agents,
                                                   n_workplaces = N_WORKPLACES,
                                                   edu_to_wp_strength = EDU_TO_WP)
  cal_beta <- c(BETA_EDU, BETA_WP)
  cal_vars <- c("edu_diff", "same_wp")
  cal <- calibrate_threshold(cal_agents, cal_beta, cal_vars,
                             target_acceptance_rate = 0.20)
  BI_THRESHOLD <- cal$threshold
  BI_ROUNDS <- min(ceiling(20 / 0.20^2), 500)
  rm(cal_agents); gc(verbose = FALSE)
  cat(sprintf("  Bilateral threshold (20%% accept): %.3f, rounds: %d\n",
              BI_THRESHOLD, BI_ROUNDS))

  for (mt in c("one_sided", "bilateral")) {
    cat(sprintf("\n\n  *** Matching type: %s ***\n", mt))

    # Set bilateral-specific params (ignored by one_sided path)
    bi_args <- list()
    if (mt == "bilateral") {
      bi_args <- list(bilateral_threshold = BI_THRESHOLD,
                      bilateral_rounds = BI_ROUNDS)
    }

    # Stage A: No confounding
    cat("\n--- Stage A: No confounding ---\n")
    if (mt == "one_sided") {
      cat("  Running one diagnostic iteration to show the DGP...\n")
      set.seed(99)
      run_module3_sim(include_confounder = FALSE, matching_type = mt, verbose = TRUE)
    }
    mc_A <- do.call(run_monte_carlo, c(list(
      sim_fn = run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = FALSE, matching_type = mt), bi_args))
    rec_A <- assess_recovery(mc_A, true_M1, "M1")
    med_A <- assess_mediation(mc_A, "edu_diff", "M0", "M1")
    cat("  M1 recovery:\n"); print(rec_A[, c("parameter", "true_value", "mean_estimate", "bias", "coverage_95")])
    cat("  Mediation:\n"); print(med_A[, c("mean_pct_reduction", "sd_pct_reduction")])

    # Stage B: Confounding (baseline)
    cat("\n--- Stage B: Confounding (baseline) ---\n")
    if (mt == "one_sided") {
      cat("  Running one diagnostic iteration to show the DGP...\n")
      set.seed(99)
      run_module3_sim(include_confounder = TRUE, matching_type = mt, verbose = TRUE)
    }
    mc_B <- do.call(run_monte_carlo, c(list(
      sim_fn = run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = TRUE, matching_type = mt), bi_args))
    rec_B <- assess_recovery(mc_B, true_M1, "M1")
    med_B <- assess_mediation(mc_B, "edu_diff", "M0", "M1")
    cat("  M1 recovery:\n"); print(rec_B[, c("parameter", "true_value", "mean_estimate", "bias", "coverage_95")])
    cat("  Mediation:\n"); print(med_B[, c("mean_pct_reduction", "sd_pct_reduction")])

    # Stage C: Confounder correlated with education
    cat("\n--- Stage C: Confounder correlated with education (rho=0.3) ---\n")
    if (mt == "one_sided") {
      cat("  Running one diagnostic iteration to show the DGP...\n")
      set.seed(99)
      run_module3_sim(include_confounder = TRUE, openness_edu_cor = 0.3,
                      matching_type = mt, verbose = TRUE)
    }
    mc_C <- do.call(run_monte_carlo, c(list(
      sim_fn = run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = TRUE, openness_edu_cor = 0.3, matching_type = mt), bi_args))
    rec_C <- assess_recovery(mc_C, true_M1, "M1")
    med_C <- assess_mediation(mc_C, "edu_diff", "M0", "M1")
    cat("  M1 recovery:\n"); print(rec_C[, c("parameter", "true_value", "mean_estimate", "bias", "coverage_95")])
    cat("  Mediation:\n"); print(med_C[, c("mean_pct_reduction", "sd_pct_reduction")])

    # Stage D: Choice set selection (one-sided only — bilateral has its own selection)
    if (mt == "one_sided") {
      cat("\n--- Stage D: Choice set selection (strength=1.0, rho=0.3) ---\n")
      cat("  Running one diagnostic iteration to show the DGP...\n")
      set.seed(99)
      run_module3_sim(include_confounder = TRUE, openness_edu_cor = 0.3,
                      openness_cs_strength = 1.0, matching_type = mt, verbose = TRUE)
      mc_D <- run_monte_carlo(
        run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
        include_confounder = TRUE, openness_edu_cor = 0.3,
        openness_cs_strength = 1.0, matching_type = mt
      )
      rec_D <- assess_recovery(mc_D, true_M1, "M1")
      med_D <- assess_mediation(mc_D, "edu_diff", "M0", "M1")
      cat("  M1 recovery:\n"); print(rec_D[, c("parameter", "true_value", "mean_estimate", "bias", "coverage_95")])
      cat("  Mediation:\n"); print(med_D[, c("mean_pct_reduction", "sd_pct_reduction")])
    } else {
      mc_D <- NULL
    }

    all_results[[mt]] <- list(A = mc_A, B = mc_B, C = mc_C, D = mc_D)
  }

  # --- Summary comparison: one-sided vs bilateral ---
  cat("\n\n--- COMPARISON: One-sided vs Bilateral (Stage B, confounding) ---\n")
  cat(sprintf("  %-25s %12s %12s %12s %12s\n",
              "", "edu_bias", "wp_bias", "med_pct", "med_sd"))
  cat("  ", rep("-", 75), "\n", sep = "")

  for (mt in c("one_sided", "bilateral")) {
    rec <- assess_recovery(all_results[[mt]]$B, true_M1, "M1")
    med <- assess_mediation(all_results[[mt]]$B, "edu_diff", "M0", "M1")
    cat(sprintf("  %-25s %12.3f %12.3f %12.1f %12.1f\n",
                mt,
                rec$bias[rec$parameter == "edu_diff"],
                rec$bias[rec$parameter == "same_wp"],
                med$mean_pct_reduction,
                med$sd_pct_reduction))
  }

  all_results
}


# =============================================================================
# PART II: ASSUMPTION SWEEPS
# =============================================================================

run_assumption_sweep <- function(assumption_name, vary_param, vary_values,
                                  fixed_params = list(), R = R_SIMS,
                                  matching_type = "one_sided") {
  #' Generic function to sweep one assumption while holding others at baseline.
  #'
  #' @param assumption_name Character label for the assumption
  #' @param vary_param Character name of the parameter to vary in run_module3_sim
  #' @param vary_values Vector of values to test
  #' @param fixed_params Named list of non-default parameters to hold fixed
  #' @param matching_type "one_sided" or "bilateral"
  #' @return data.frame with recovery and mediation stats per level

  cat(sprintf("\n--- Assumption sweep: %s [%s] ---\n", assumption_name, matching_type))

  sweep_results <- list()

  for (val in vary_values) {
    cat(sprintf("  %s = %s\n", vary_param, as.character(val)))

    # Build parameter list
    params <- c(fixed_params, setNames(list(val), vary_param))
    params$matching_type <- matching_type
    params$R <- R
    params$n_cores <- N_CORES
    params$seed <- 42
    params$sim_fn <- run_module3_sim

    mc <- do.call(run_monte_carlo, params)

    # Update true values if we're sweeping a DGP coefficient
    true_M1 <- c(edu_diff = BETA_EDU, same_wp = BETA_WP)
    if (vary_param == "beta_edu") true_M1["edu_diff"] <- val
    if (vary_param == "beta_wp") true_M1["same_wp"] <- val

    rec <- assess_recovery(mc, true_M1, "M1")
    med <- assess_mediation(mc, "edu_diff", "M0", "M1")

    row <- data.frame(
      assumption = assumption_name,
      matching_type = matching_type,
      level = as.character(val),
      level_numeric = as.numeric(val),
      edu_diff_bias = rec$bias[rec$parameter == "edu_diff"],
      edu_diff_rmse = rec$rmse[rec$parameter == "edu_diff"],
      edu_diff_coverage = rec$coverage_95[rec$parameter == "edu_diff"],
      same_wp_bias = rec$bias[rec$parameter == "same_wp"],
      same_wp_rmse = rec$rmse[rec$parameter == "same_wp"],
      same_wp_coverage = rec$coverage_95[rec$parameter == "same_wp"],
      mediation_pct = med$mean_pct_reduction,
      mediation_sd = med$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    sweep_results[[as.character(val)]] <- row
    cat(sprintf("    edu bias=%.3f, wp bias=%.3f, mediation=%.1f%%\n",
                row$edu_diff_bias, row$same_wp_bias, row$mediation_pct))
  }

  bind_rows(sweep_results)
}


run_part_II <- function() {
  cat("\n")
  cat("=" , rep("=", 70), "\n", sep = "")
  cat("PART II: Assumption Sweeps\n")
  cat("=" , rep("=", 70), "\n")

  all_sweeps <- list()

  for (mt in c("one_sided")) {  # Oracle only — bilateral tested in Module 4
    cat(sprintf("\n\n  *** Matching type: %s ***\n", mt))
    mt_prefix <- if (mt == "bilateral") "bi_" else ""

    # --- A1: Mediation pathway strength ---
    all_sweeps[[paste0(mt_prefix, "A1")]] <- run_assumption_sweep(
      "A1: Mediation pathway strength",
      vary_param = "edu_to_wp",
      vary_values = c(0.5, 1.0, 2.0, 3.0, 5.0),
      fixed_params = list(include_confounder = TRUE),
      matching_type = mt
    )

    # --- A2: Choice set size ---
    all_sweeps[[paste0(mt_prefix, "A2")]] <- run_assumption_sweep(
      "A2: Choice set size (J)",
      vary_param = "J",
      vary_values = c(10, 20, 30, 50, 100),
      fixed_params = list(include_confounder = TRUE),
      matching_type = mt
    )

    # --- A3: Confounder type ---
    cat(sprintf("\n--- Assumption sweep: A3: Confounder type [%s] ---\n", mt))
    true_M1 <- c(edu_diff = BETA_EDU, same_wp = BETA_WP)
    A3_results <- list()

    for (ctype in c("homophily", "level", "interaction")) {
      cat(sprintf("  confounder_type = %s\n", ctype))
      mc <- run_monte_carlo(
        run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
        include_confounder = TRUE, confounder_type = ctype,
        matching_type = mt
      )
      rec <- assess_recovery(mc, true_M1, "M1")
      med <- assess_mediation(mc, "edu_diff", "M0", "M1")

      A3_results[[ctype]] <- data.frame(
        assumption = "A3: Confounder type",
        matching_type = mt,
        level = ctype,
        level_numeric = match(ctype, c("homophily", "level", "interaction")),
        edu_diff_bias = rec$bias[rec$parameter == "edu_diff"],
        edu_diff_rmse = rec$rmse[rec$parameter == "edu_diff"],
        edu_diff_coverage = rec$coverage_95[rec$parameter == "edu_diff"],
        same_wp_bias = rec$bias[rec$parameter == "same_wp"],
        same_wp_rmse = rec$rmse[rec$parameter == "same_wp"],
        same_wp_coverage = rec$coverage_95[rec$parameter == "same_wp"],
        mediation_pct = med$mean_pct_reduction,
        mediation_sd = med$sd_pct_reduction,
        stringsAsFactors = FALSE
      )
      cat(sprintf("    edu bias=%.3f, wp bias=%.3f, mediation=%.1f%%\n",
                  A3_results[[ctype]]$edu_diff_bias,
                  A3_results[[ctype]]$same_wp_bias,
                  A3_results[[ctype]]$mediation_pct))
    }
    all_sweeps[[paste0(mt_prefix, "A3")]] <- bind_rows(A3_results)

    # --- A4: Group structure ---
    cat(sprintf("\n--- Assumption sweep: A4: Group structure [%s] ---\n", mt))
    A4_results <- list()

    # No groups (baseline)
    cat("  No groups (baseline)\n")
    mc_nogrp <- run_monte_carlo(
      run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = TRUE, use_groups = FALSE, matching_type = mt
    )
    rec_ng <- assess_recovery(mc_nogrp, true_M1, "M1")
    med_ng <- assess_mediation(mc_nogrp, "edu_diff", "M0", "M1")
    A4_results[["no_groups"]] <- data.frame(
      assumption = "A4: Group structure", matching_type = mt,
      level = "no_groups", level_numeric = 0,
      edu_diff_bias = rec_ng$bias[rec_ng$parameter == "edu_diff"],
      edu_diff_rmse = rec_ng$rmse[rec_ng$parameter == "edu_diff"],
      edu_diff_coverage = rec_ng$coverage_95[rec_ng$parameter == "edu_diff"],
      same_wp_bias = rec_ng$bias[rec_ng$parameter == "same_wp"],
      same_wp_rmse = rec_ng$rmse[rec_ng$parameter == "same_wp"],
      same_wp_coverage = rec_ng$coverage_95[rec_ng$parameter == "same_wp"],
      mediation_pct = med_ng$mean_pct_reduction,
      mediation_sd = med_ng$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    # With groups, no group controls in model
    cat("  With groups (70/20/10), no group controls\n")
    mc_grp <- run_monte_carlo(
      run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = TRUE, use_groups = TRUE, matching_type = mt
    )
    rec_g <- assess_recovery(mc_grp, true_M1, "M1")
    med_g <- assess_mediation(mc_grp, "edu_diff", "M0", "M1")
    A4_results[["groups_no_ctrl"]] <- data.frame(
      assumption = "A4: Group structure", matching_type = mt,
      level = "groups_no_ctrl", level_numeric = 1,
      edu_diff_bias = rec_g$bias[rec_g$parameter == "edu_diff"],
      edu_diff_rmse = rec_g$rmse[rec_g$parameter == "edu_diff"],
      edu_diff_coverage = rec_g$coverage_95[rec_g$parameter == "edu_diff"],
      same_wp_bias = rec_g$bias[rec_g$parameter == "same_wp"],
      same_wp_rmse = rec_g$rmse[rec_g$parameter == "same_wp"],
      same_wp_coverage = rec_g$coverage_95[rec_g$parameter == "same_wp"],
      mediation_pct = med_g$mean_pct_reduction,
      mediation_sd = med_g$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    # With groups + group controls
    cat("  With groups (70/20/10), with group controls\n")
    true_M1g <- c(edu_diff = BETA_EDU, same_wp = BETA_WP)
    rec_gc <- assess_recovery(mc_grp, true_M1g, "M1g")
    med_gc <- assess_mediation(mc_grp, "edu_diff", "M0g", "M1g")
    A4_results[["groups_with_ctrl"]] <- data.frame(
      assumption = "A4: Group structure", matching_type = mt,
      level = "groups_with_ctrl", level_numeric = 2,
      edu_diff_bias = rec_gc$bias[rec_gc$parameter == "edu_diff"],
      edu_diff_rmse = rec_gc$rmse[rec_gc$parameter == "edu_diff"],
      edu_diff_coverage = rec_gc$coverage_95[rec_gc$parameter == "edu_diff"],
      same_wp_bias = rec_gc$bias[rec_gc$parameter == "same_wp"],
      same_wp_rmse = rec_gc$rmse[rec_gc$parameter == "same_wp"],
      same_wp_coverage = rec_gc$coverage_95[rec_gc$parameter == "same_wp"],
      mediation_pct = med_gc$mean_pct_reduction,
      mediation_sd = med_gc$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    for (nm in names(A4_results)) {
      cat(sprintf("    %s: edu bias=%.3f, wp bias=%.3f, mediation=%.1f%%\n",
                  nm, A4_results[[nm]]$edu_diff_bias,
                  A4_results[[nm]]$same_wp_bias,
                  A4_results[[nm]]$mediation_pct))
    }
    all_sweeps[[paste0(mt_prefix, "A4")]] <- bind_rows(A4_results)

    # --- A5: Multiple mediators ---
    cat(sprintf("\n--- Assumption sweep: A5: Multiple mediators [%s] ---\n", mt))
    A5_results <- list()

    # Single mediator (baseline) — reuse from A4
    A5_results[["single"]] <- data.frame(
      assumption = "A5: Multiple mediators", matching_type = mt,
      level = "single", level_numeric = 1,
      edu_diff_bias = rec_ng$bias[rec_ng$parameter == "edu_diff"],
      edu_diff_rmse = rec_ng$rmse[rec_ng$parameter == "edu_diff"],
      edu_diff_coverage = rec_ng$coverage_95[rec_ng$parameter == "edu_diff"],
      same_wp_bias = rec_ng$bias[rec_ng$parameter == "same_wp"],
      same_wp_rmse = rec_ng$rmse[rec_ng$parameter == "same_wp"],
      same_wp_coverage = rec_ng$coverage_95[rec_ng$parameter == "same_wp"],
      mediation_pct = med_ng$mean_pct_reduction,
      mediation_sd = med_ng$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    # Two mediators, confounder affects only workplace
    cat("  Two mediators, confounder -> workplace only\n")
    mc_2med_wp <- run_monte_carlo(
      run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = TRUE, use_neighbourhood = TRUE,
      openness_to_nbhd = 0.0, matching_type = mt
    )
    true_M2 <- c(edu_diff = BETA_EDU, same_wp = BETA_WP, same_nbhd = 1.0)
    rec_2a <- assess_recovery(mc_2med_wp, true_M2, "M2")
    med_2a_total <- assess_mediation(mc_2med_wp, "edu_diff", "M0", "M2")
    A5_results[["two_med_wp_conf"]] <- data.frame(
      assumption = "A5: Multiple mediators", matching_type = mt,
      level = "two_wp_conf", level_numeric = 2,
      edu_diff_bias = rec_2a$bias[rec_2a$parameter == "edu_diff"],
      edu_diff_rmse = rec_2a$rmse[rec_2a$parameter == "edu_diff"],
      edu_diff_coverage = rec_2a$coverage_95[rec_2a$parameter == "edu_diff"],
      same_wp_bias = rec_2a$bias[rec_2a$parameter == "same_wp"],
      same_wp_rmse = rec_2a$rmse[rec_2a$parameter == "same_wp"],
      same_wp_coverage = rec_2a$coverage_95[rec_2a$parameter == "same_wp"],
      mediation_pct = med_2a_total$mean_pct_reduction,
      mediation_sd = med_2a_total$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    # Two mediators, confounder affects both
    cat("  Two mediators, confounder -> both\n")
    mc_2med_both <- run_monte_carlo(
      run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
      include_confounder = TRUE, use_neighbourhood = TRUE,
      openness_to_nbhd = 1.0, matching_type = mt
    )
    rec_2b <- assess_recovery(mc_2med_both, true_M2, "M2")
    med_2b_total <- assess_mediation(mc_2med_both, "edu_diff", "M0", "M2")
    A5_results[["two_med_both_conf"]] <- data.frame(
      assumption = "A5: Multiple mediators", matching_type = mt,
      level = "two_both_conf", level_numeric = 3,
      edu_diff_bias = rec_2b$bias[rec_2b$parameter == "edu_diff"],
      edu_diff_rmse = rec_2b$rmse[rec_2b$parameter == "edu_diff"],
      edu_diff_coverage = rec_2b$coverage_95[rec_2b$parameter == "edu_diff"],
      same_wp_bias = rec_2b$bias[rec_2b$parameter == "same_wp"],
      same_wp_rmse = rec_2b$rmse[rec_2b$parameter == "same_wp"],
      same_wp_coverage = rec_2b$coverage_95[rec_2b$parameter == "same_wp"],
      mediation_pct = med_2b_total$mean_pct_reduction,
      mediation_sd = med_2b_total$sd_pct_reduction,
      stringsAsFactors = FALSE
    )

    for (nm in names(A5_results)) {
      cat(sprintf("    %s: edu bias=%.3f, wp bias=%.3f, mediation=%.1f%%\n",
                  nm, A5_results[[nm]]$edu_diff_bias,
                  A5_results[[nm]]$same_wp_bias,
                  A5_results[[nm]]$mediation_pct))
    }
    all_sweeps[[paste0(mt_prefix, "A5")]] <- bind_rows(A5_results)

    # --- A6: Workplace concentration ---
    cat(sprintf("\n--- Assumption sweep: A6: Workplace concentration [%s] ---\n", mt))
    A6_results <- list()

    for (conc in c("equal", "zipf")) {
      cat(sprintf("  concentration = %s\n", conc))
      mc <- run_monte_carlo(
        run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
        include_confounder = TRUE, wp_concentration = conc, matching_type = mt
      )
      rec <- assess_recovery(mc, true_M1, "M1")
      med <- assess_mediation(mc, "edu_diff", "M0", "M1")
      A6_results[[conc]] <- data.frame(
        assumption = "A6: Workplace concentration", matching_type = mt,
        level = conc, level_numeric = match(conc, c("equal", "zipf")),
        edu_diff_bias = rec$bias[rec$parameter == "edu_diff"],
        edu_diff_rmse = rec$rmse[rec$parameter == "edu_diff"],
        edu_diff_coverage = rec$coverage_95[rec$parameter == "edu_diff"],
        same_wp_bias = rec$bias[rec$parameter == "same_wp"],
        same_wp_rmse = rec$rmse[rec$parameter == "same_wp"],
        same_wp_coverage = rec$coverage_95[rec$parameter == "same_wp"],
        mediation_pct = med$mean_pct_reduction,
        mediation_sd = med$sd_pct_reduction,
        stringsAsFactors = FALSE
      )
      cat(sprintf("    edu bias=%.3f, wp bias=%.3f, mediation=%.1f%%\n",
                  A6_results[[conc]]$edu_diff_bias,
                  A6_results[[conc]]$same_wp_bias,
                  A6_results[[conc]]$mediation_pct))
    }
    all_sweeps[[paste0(mt_prefix, "A6")]] <- bind_rows(A6_results)

    # --- A7: Direct vs indirect effect ratio ---
    all_sweeps[[paste0(mt_prefix, "A7")]] <- run_assumption_sweep(
      "A7: Direct effect magnitude (beta_edu)",
      vary_param = "beta_edu",
      vary_values = c(-0.3, -0.5, -1.0, -2.0, -3.0),
      fixed_params = list(include_confounder = TRUE),
      matching_type = mt
    )

    # --- A8: Confounding strength sweep ---
    all_sweeps[[paste0(mt_prefix, "A8")]] <- run_assumption_sweep(
      "A8: Confounding strength (openness -> wp)",
      vary_param = "openness_to_wp",
      vary_values = c(0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0),
      fixed_params = list(include_confounder = TRUE),
      matching_type = mt
    )

  }  # end matching_type loop

  return(all_sweeps)
}


# =============================================================================
# PART III: SELECTED INTERACTIONS
# =============================================================================

run_part_III <- function() {
  cat("\n")
  cat("=" , rep("=", 70), "\n", sep = "")
  cat("PART III: Interaction Tests\n")
  cat("=" , rep("=", 70), "\n")

  true_M1 <- c(edu_diff = BETA_EDU, same_wp = BETA_WP)

  all_interactions <- list()

  for (mt in c("one_sided")) {  # Oracle only — bilateral tested in Module 4
    cat(sprintf("\n\n  *** Matching type: %s ***\n", mt))

    # --- Interaction 1: Confounding strength × Choice set size ---
    cat("\n--- Confounding strength x Choice set size ---\n")
    int1_results <- list()

    for (cs in c(0.5, 1.0, 2.0)) {
      for (J in c(10, 30, 100)) {
        label <- sprintf("conf=%.1f_J=%d", cs, J)
        cat(sprintf("  %s\n", label))

        mc <- run_monte_carlo(
          run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
          include_confounder = TRUE, openness_to_wp = cs, J = J,
          matching_type = mt
        )
        rec <- assess_recovery(mc, true_M1, "M1")
        med <- assess_mediation(mc, "edu_diff", "M0", "M1")

        int1_results[[label]] <- data.frame(
          matching_type = mt,
          confounding = cs, J = J,
          edu_diff_bias = rec$bias[rec$parameter == "edu_diff"],
          same_wp_bias = rec$bias[rec$parameter == "same_wp"],
          mediation_pct = med$mean_pct_reduction,
          stringsAsFactors = FALSE
        )
      }
    }

    int1_df <- bind_rows(int1_results)
    cat("\n  Results:\n")
    print(int1_df)

    # --- Interaction 2: Confounding strength × Group structure ---
    cat("\n--- Confounding strength x Group structure ---\n")
    int2_results <- list()

    for (cs in c(0.5, 1.0, 2.0)) {
      for (grp in c(FALSE, TRUE)) {
        label <- sprintf("conf=%.1f_groups=%s", cs, grp)
        cat(sprintf("  %s\n", label))

        mc <- run_monte_carlo(
          run_module3_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
          include_confounder = TRUE, openness_to_wp = cs, use_groups = grp,
          matching_type = mt
        )
        rec <- assess_recovery(mc, true_M1, "M1")
        med <- assess_mediation(mc, "edu_diff", "M0", "M1")

        int2_results[[label]] <- data.frame(
          matching_type = mt,
          confounding = cs, groups = grp,
          edu_diff_bias = rec$bias[rec$parameter == "edu_diff"],
          same_wp_bias = rec$bias[rec$parameter == "same_wp"],
          mediation_pct = med$mean_pct_reduction,
          stringsAsFactors = FALSE
        )
      }
    }

    int2_df <- bind_rows(int2_results)
    cat("\n  Results:\n")
    print(int2_df)

    mt_key <- if (mt == "bilateral") "bi_" else ""
    all_interactions[[paste0(mt_key, "conf_x_J")]] <- int1_df
    all_interactions[[paste0(mt_key, "conf_x_groups")]] <- int2_df
  }

  all_interactions
}


# =============================================================================
# PLOTS
# =============================================================================

generate_plots <- function(sweeps, interactions) {
  cat("\n--- Generating plots ---\n")
  dir.create("results", showWarnings = FALSE)

  all_sweeps <- bind_rows(sweeps)

  # --- Master sweep plot: bias in beta_wp by assumption and matching type ---
  if ("matching_type" %in% names(all_sweeps)) {
    p_wp_bias <- ggplot(all_sweeps, aes(x = level, y = same_wp_bias,
                                         fill = matching_type)) +
      geom_col(position = "dodge", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ assumption, scales = "free_x", ncol = 2) +
      labs(title = "Bias in Workplace Coefficient Across Assumptions",
           x = "", y = "Bias (beta_wp estimate - true)",
           fill = "Matching") +
      scale_fill_manual(values = c("one_sided" = "steelblue",
                                    "bilateral" = "coral")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  } else {
    p_wp_bias <- ggplot(all_sweeps, aes(x = level, y = same_wp_bias)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ assumption, scales = "free_x", ncol = 2) +
      labs(title = "Bias in Workplace Coefficient Across Assumptions",
           x = "", y = "Bias (beta_wp estimate - true)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  }
  ggsave("results/module3_wp_bias_all_assumptions.pdf", p_wp_bias, width = 14, height = 12)

  # --- Master sweep plot: mediation share ---
  if ("matching_type" %in% names(all_sweeps)) {
    p_mediation <- ggplot(all_sweeps, aes(x = level, y = mediation_pct,
                                           fill = matching_type)) +
      geom_col(position = "dodge", alpha = 0.7) +
      facet_wrap(~ assumption, scales = "free_x", ncol = 2) +
      labs(title = "Estimated Mediation Share Across Assumptions",
           x = "", y = "Mediation share (%)", fill = "Matching") +
      scale_fill_manual(values = c("one_sided" = "darkorange",
                                    "bilateral" = "purple")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  } else {
    p_mediation <- ggplot(all_sweeps, aes(x = level, y = mediation_pct)) +
      geom_col(fill = "darkorange", alpha = 0.7) +
      facet_wrap(~ assumption, scales = "free_x", ncol = 2) +
      labs(title = "Estimated Mediation Share Across Assumptions",
           x = "", y = "Mediation share (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  }
  ggsave("results/module3_mediation_all_assumptions.pdf", p_mediation, width = 14, height = 12)

  # --- Confounding strength (A8) detailed: one_sided vs bilateral ---
  a8_keys <- grep("A8$", names(sweeps), value = TRUE)
  if (length(a8_keys) > 0) {
    a8 <- bind_rows(sweeps[a8_keys])
    if ("matching_type" %in% names(a8)) {
      p_a8 <- ggplot(a8, aes(x = level_numeric, linetype = matching_type)) +
        geom_line(aes(y = same_wp_bias, color = "beta_wp bias"), linewidth = 1) +
        geom_line(aes(y = edu_diff_bias, color = "beta_edu bias"), linewidth = 1) +
        geom_point(aes(y = same_wp_bias, color = "beta_wp bias"), size = 3) +
        geom_point(aes(y = edu_diff_bias, color = "beta_edu bias"), size = 3) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(title = "A8: Parameter Bias by Confounding Strength",
             x = "Openness -> Workplace strength", y = "Bias",
             color = "", linetype = "Matching") +
        theme_minimal()
    } else {
      p_a8 <- ggplot(a8, aes(x = level_numeric)) +
        geom_line(aes(y = same_wp_bias, color = "beta_wp bias"), linewidth = 1) +
        geom_line(aes(y = edu_diff_bias, color = "beta_edu bias"), linewidth = 1) +
        geom_point(aes(y = same_wp_bias, color = "beta_wp bias"), size = 3) +
        geom_point(aes(y = edu_diff_bias, color = "beta_edu bias"), size = 3) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(title = "A8: Parameter Bias by Confounding Strength",
             x = "Openness -> Workplace strength", y = "Bias", color = "") +
        theme_minimal()
    }
    ggsave("results/module3_A8_confounding_bias.pdf", p_a8, width = 8, height = 5)
  }

  # --- Interaction plots: combine matching types ---
  conf_J_keys <- grep("conf_x_J$", names(interactions), value = TRUE)
  if (length(conf_J_keys) > 0) {
    int1_df <- bind_rows(interactions[conf_J_keys])
    p_int1 <- ggplot(int1_df,
                     aes(x = factor(J), y = same_wp_bias, fill = factor(confounding))) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ matching_type) +
      labs(title = "Interaction: Confounding x Choice Set Size",
           x = "Choice set size (J)", y = "beta_wp bias",
           fill = "Confounding\nstrength") +
      theme_minimal()
    ggsave("results/module3_interaction_conf_J.pdf", p_int1, width = 10, height = 5)
  }

  conf_grp_keys <- grep("conf_x_groups$", names(interactions), value = TRUE)
  if (length(conf_grp_keys) > 0) {
    int2_df <- bind_rows(interactions[conf_grp_keys])
    p_int2 <- ggplot(int2_df,
                     aes(x = factor(confounding), y = same_wp_bias, fill = groups)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ matching_type) +
      labs(title = "Interaction: Confounding x Group Structure",
           x = "Confounding strength", y = "beta_wp bias",
           fill = "Groups") +
      theme_minimal()
    ggsave("results/module3_interaction_conf_groups.pdf", p_int2, width = 10, height = 5)
  }

  cat("Plots saved to results/\n")
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("=" , rep("=", 70), "\n", sep = "")
cat("MODULE 3: Mediation — Expanded Assumption Testing\n")
cat("  (with one-sided and bilateral matching comparison)\n")
cat("=" , rep("=", 70), "\n")
cat(sprintf("Settings: N=%d, J=%d (baseline), R=%d reps per condition\n",
            N_AGENTS, J_ALTS, R_SIMS))

dir.create("results", showWarnings = FALSE)

# Run all parts (each iterates over one_sided and bilateral)
part_I <- run_part_I()
part_II <- run_part_II()
part_III <- run_part_III()

# Generate plots
generate_plots(part_II, part_III)

# =============================================================================
# SAVE CSV RESULTS
# =============================================================================

# --- Part I: Core stages A-D ---
true_M1 <- c(edu_diff = BETA_EDU, same_wp = BETA_WP)
part_I_rows <- list()
for (mt in names(part_I)) {
  for (stage in names(part_I[[mt]])) {
    mc <- part_I[[mt]][[stage]]
    if (is.null(mc)) next
    rec <- assess_recovery(mc, true_M1, "M1")
    med <- assess_mediation(mc, "edu_diff", "M0", "M1")
    for (i in seq_len(nrow(rec))) {
      part_I_rows[[length(part_I_rows) + 1]] <- data.frame(
        part = "I",
        matching_type = mt,
        stage = stage,
        parameter = rec$parameter[i],
        true_value = rec$true_value[i],
        mean_est = rec$mean_estimate[i],
        bias = rec$bias[i],
        rmse = rec$rmse[i],
        coverage_95 = rec$coverage_95[i],
        mediation_pct = med$mean_pct_reduction,
        mediation_sd = med$sd_pct_reduction,
        stringsAsFactors = FALSE
      )
    }
  }
}
part_I_df <- bind_rows(part_I_rows)

# --- Part II: Assumption sweeps ---
part_II_df <- bind_rows(part_II)
part_II_df$part <- "II"

# --- Part III: Interactions ---
part_III_df <- bind_rows(part_III)
part_III_df$part <- "III"

# Save each part
write.csv(part_I_df, "results/module3_core_stages.csv", row.names = FALSE)
write.csv(part_II_df, "results/module3_assumption_sweeps.csv", row.names = FALSE)
write.csv(part_III_df, "results/module3_interactions.csv", row.names = FALSE)

cat("\nCSV results saved to results/module3_*.csv\n")
cat("Plots saved to results/\n")
cat("--- Module 3 complete ---\n")
