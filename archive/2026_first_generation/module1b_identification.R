# =============================================================================
# Module 1b: Identification Tests
# =============================================================================
#
# Two sets of questions:
#
# PART A — HOMOPHILY VS STATUS MAXIMISATION
#   Can clogit distinguish between two different generating mechanisms that
#   produce similar patterns in the data?
#
#   Homophily: people prefer partners similar to themselves on education.
#     U(i,j) = beta_h * |edu_i - edu_j| + beta_sg * same_group + eps
#
#   Status maximisation: everyone prefers the highest-educated partner they
#     can get. The resulting "winnowing" produces assortative matching that
#     looks like homophily but is driven by a hierarchy, not similarity.
#     U(i,j) = beta_s * alt_edu + beta_sg * same_group + eps
#
#   Mixed: both mechanisms operate simultaneously.
#     U(i,j) = beta_h * |edu_diff| + beta_s * alt_edu + beta_sg * same_group + eps
#
#   We generate data under each DGP and estimate multiple clogit specifications
#   to see which can correctly identify the true mechanism.
#
# PART B — CONDITIONAL LOGIT VS MULTINOMIAL LOGIT
#   When is the multinomial logit (mlogit) equivalent to clogit for studying
#   partner choice, and when does clogit offer an advantage?
#
#   The mlogit models P(partner_group = j | ego characteristics). It treats
#   partner groups as discrete categories and does not condition on individual
#   choice sets. The clogit conditions on individual-level choice sets.
#
#   We test:
#   B1: Baseline equivalence (balanced groups, no opportunity variation)
#   B2: Group size sensitivity (do unequal groups bias the mlogit?)
#   B3: Opportunity variation (local markets with different group composition)
#   B4: Can mlogit be corrected with composition controls?
#
# =============================================================================

source("utils.R")
library(nnet)  # for multinomial logit


# =============================================================================
# GLOBAL PARAMETERS
# =============================================================================

R_SIMS <- 200
N_AGENTS <- 1000
J_ALTS <- 30
N_CORES <- 1

# Part A DGP parameters
BETA_HOMOPHILY <- -1.5    # homophily: prefer similar education
BETA_STATUS <- 0.5        # status max: prefer high-educated partner
BETA_SAME_GROUP <- 0.8    # endogamy: prefer same group

# Group structure (used throughout)
GROUP_PROPS <- c(0.70, 0.20, 0.10)
EDU_MEANS <- c(12, 10, 14)
EDU_SDS <- c(2, 3, 1.5)


# =============================================================================
# PART A: HOMOPHILY VS STATUS MAXIMISATION
# =============================================================================

# =============================================================================
# A: SIMULATION FUNCTION
# =============================================================================

run_identification_sim <- function(
    N = N_AGENTS,
    J = J_ALTS,
    group_proportions = GROUP_PROPS,
    edu_means = EDU_MEANS,
    edu_sds = EDU_SDS,
    # DGP type: "homophily", "status_max", "mixed"
    dgp_type = "homophily",
    # Matching type: "one_sided" (original) or "bilateral"
    matching_type = "one_sided",
    # Bilateral matching parameters
    bilateral_threshold = -0.5,
    bilateral_rounds = 20,
    beta_h = BETA_HOMOPHILY,
    beta_s = BETA_STATUS,
    beta_sg = BETA_SAME_GROUP,
    verbose = FALSE
) {
  n_groups <- length(group_proportions)

  # 1. Generate agents
  agents <- generate_agents(
    N = N, n_groups = n_groups,
    group_proportions = group_proportions,
    edu_means = edu_means, edu_sds = edu_sds
  )

  # Define DGP parameters
  if (dgp_type == "homophily") {
    gen_beta <- c(beta_h, beta_sg)
    gen_vars <- c("edu_diff", "same_group")
  } else if (dgp_type == "status_max") {
    gen_beta <- c(beta_s, beta_sg)
    gen_vars <- c("alt_edu", "same_group")
  } else if (dgp_type == "mixed") {
    gen_beta <- c(beta_h, beta_s, beta_sg)
    gen_vars <- c("edu_diff", "alt_edu", "same_group")
  }

  # 2. Generate unions and construct choice sets
  if (matching_type == "bilateral") {
    # --- Two-sided matching ---
    bm <- bilateral_matching(
      agents, beta = gen_beta, var_names = gen_vars,
      n_rounds = bilateral_rounds, threshold = bilateral_threshold,
      verbose = FALSE
    )

    if (nrow(bm$matches) < 20) {
      warning("Bilateral matching produced too few pairs")
      return(NULL)
    }

    # Construct choice sets from matched pairs
    cd <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    cd <- compute_match_variables(cd)

    # Add alt_edu and interaction variables
    cd$alt_edu <- cd$alt_education
    cd$ego_x_alt_edu <- cd$ego_education * cd$alt_education

    # Diagnostics for bilateral
    if (verbose) {
      cat(sprintf("\n  DGP: %s (BILATERAL MATCHING)\n", dgp_type))
      cat(sprintf("  Pairing rate: %.1f%% (%d pairs, %d unmatched)\n",
                  100 * bm$pairing_rate, nrow(bm$matches), length(bm$unmatched)))
      describe_population(agents)

      # Union patterns from actual matches
      m <- bm$matches
      ego_edu <- agents$education[match(m$ego_id, agents$id)]
      alt_edu_v <- agents$education[match(m$alt_id, agents$id)]
      ego_grp <- agents$group[match(m$ego_id, agents$id)]
      alt_grp <- agents$group[match(m$alt_id, agents$id)]

      cat("\n  ---- BILATERAL UNION PATTERNS ----\n")
      cat(sprintf("  Same group rate: %.3f\n", mean(ego_grp == alt_grp)))
      cat(sprintf("  Mean |edu_diff|: %.2f\n", mean(abs(ego_edu - alt_edu_v))))
      cat(sprintf("  Mean alt_edu: %.2f\n", mean(alt_edu_v)))
      cat(sprintf("  Cor(ego_edu, alt_edu): %.3f\n", cor(ego_edu, alt_edu_v)))
    }

  } else {
    # --- One-sided matching (original) ---
    cd <- construct_choice_sets(agents, J = J)
    cd <- compute_match_variables(cd)

    # Add alt_edu and interaction variables
    cd$alt_edu <- cd$alt_education
    cd$ego_x_alt_edu <- cd$ego_education * cd$alt_education

    # Simulate choices
    cd <- simulate_choices(cd, beta = gen_beta, var_names = gen_vars)

    # Diagnostics for one-sided
    if (verbose) {
      cat(sprintf("\n  DGP: %s (ONE-SIDED MATCHING)\n", dgp_type))
      describe_population(agents)

      chosen <- cd[cd$chosen == 1, ]
      cat("\n  ---- UNION PATTERNS ----\n")
      cat(sprintf("  Same group rate: %.3f\n", mean(chosen$same_group)))
      cat(sprintf("  Mean |edu_diff| among chosen: %.2f (baseline: %.2f)\n",
                  mean(chosen$edu_diff), mean(cd$edu_diff)))
      cat(sprintf("  Mean alt_edu among chosen: %.2f (baseline: %.2f)\n",
                  mean(chosen$alt_edu), mean(cd$alt_edu)))
      cat(sprintf("  Cor(ego_edu, alt_edu) among pairs: %.3f\n",
                  cor(chosen$ego_education, chosen$alt_education)))
    }
  }

  # 3. Estimate four clogit specifications
  # Spec H: homophily specification
  spec_H <- estimate_clogit(cd, "edu_diff + same_group")

  # Spec S: status maximisation specification
  spec_S <- estimate_clogit(cd, "alt_edu + same_group")

  # Spec M: mixed specification (both edu_diff and alt_edu)
  spec_M <- estimate_clogit(cd, "edu_diff + alt_edu + same_group")

  # Spec I: interaction specification (ego_edu * alt_edu captures assortative matching)
  spec_I <- estimate_clogit(cd, "ego_x_alt_edu + same_group")

  list(
    spec_H = spec_H,
    spec_S = spec_S,
    spec_M = spec_M,
    spec_I = spec_I
  )
}


# =============================================================================
# A: RUN IDENTIFICATION TESTS
# =============================================================================

get_true_values <- function(param_names, dgp) {
  #' Helper: build named vector of true parameter values given DGP and spec.
  true_vals <- numeric(length(param_names))
  names(true_vals) <- param_names
  for (p in param_names) {
    if (p == "edu_diff" && dgp %in% c("homophily", "mixed")) {
      true_vals[p] <- BETA_HOMOPHILY
    } else if (p == "alt_edu" && dgp %in% c("status_max", "mixed")) {
      true_vals[p] <- BETA_STATUS
    } else if (p == "same_group") {
      true_vals[p] <- BETA_SAME_GROUP
    } else if (p == "edu_diff" && dgp == "status_max") {
      true_vals[p] <- 0
    } else if (p == "alt_edu" && dgp == "homophily") {
      true_vals[p] <- 0
    } else if (p == "ego_x_alt_edu") {
      true_vals[p] <- 0
    } else {
      true_vals[p] <- 0
    }
  }
  return(true_vals)
}


run_part_A <- function() {
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("PART A: Homophily vs Status Maximisation — Identification\n")
  cat("=", rep("=", 70), "\n")

  dgp_types <- c("homophily", "status_max", "mixed")
  # Use oracle (one-sided) only — bilateral matching tested in Module 4
  matching_types <- c("one_sided")
  spec_names <- c("spec_H", "spec_S", "spec_M", "spec_I")

  all_A <- list()

  for (mtype in matching_types) {
    cat(sprintf("\n\n========== MATCHING TYPE: %s ==========\n",
                toupper(mtype)))

    for (dgp in dgp_types) {
      label <- paste(dgp, mtype, sep = "_")
      cat(sprintf("\n\n--- DGP: %s (%s) ---\n", toupper(dgp), mtype))

      # Diagnostic run
      cat("  Running one diagnostic iteration...\n")
      set.seed(99)
      run_identification_sim(dgp_type = dgp, matching_type = mtype,
                             verbose = TRUE)

      # MC runs
      mc <- run_monte_carlo(
        run_identification_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
        dgp_type = dgp, matching_type = mtype
      )

      # Extract results for each specification
      cat(sprintf("\n  --- Results: DGP = %s (%s) ---\n", dgp, mtype))

      for (spec in spec_names) {
        first_valid <- mc[[which(!sapply(mc, function(x) is.null(x[[spec]])))[1]]]
        if (is.null(first_valid)) next

        param_names <- names(first_valid[[spec]]$coefficients)
        if (any(is.na(param_names))) next

        true_vals <- get_true_values(param_names, dgp)

        rec <- assess_recovery(mc, true_vals, spec)
        cat(sprintf("\n  %s:\n", spec))
        print(rec[, c("parameter", "true_value", "mean_estimate", "bias",
                       "rel_bias_pct", "coverage_95")])
      }

      all_A[[label]] <- mc
    }
  }

  # Summary comparison: one-sided vs bilateral
  cat("\n\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("SUMMARY: One-sided vs Bilateral Matching — Key Differences\n")
  cat("=", rep("=", 70), "\n\n")
  cat("For each DGP, compare the CORRECTLY SPECIFIED model:\n")
  cat("  homophily -> spec_H, status_max -> spec_S, mixed -> spec_M\n\n")

  for (dgp in dgp_types) {
    correct_spec <- switch(dgp,
      homophily = "spec_H", status_max = "spec_S", mixed = "spec_M"
    )

    cat(sprintf("--- %s (spec: %s) ---\n", toupper(dgp), correct_spec))

    for (mtype in matching_types) {
      label <- paste(dgp, mtype, sep = "_")
      mc <- all_A[[label]]
      if (is.null(mc)) next

      first_valid <- mc[[which(!sapply(mc, function(x)
        is.null(x[[correct_spec]])))[1]]]
      if (is.null(first_valid)) next

      param_names <- names(first_valid[[correct_spec]]$coefficients)
      true_vals <- get_true_values(param_names, dgp)
      rec <- assess_recovery(mc, true_vals, correct_spec)

      cat(sprintf("  %s: ", mtype))
      for (i in 1:nrow(rec)) {
        cat(sprintf("%s: est=%.3f bias=%.3f  ",
                    rec$parameter[i], rec$mean_estimate[i], rec$bias[i]))
      }
      cat("\n")
    }
    cat("\n")
  }

  return(all_A)
}


# =============================================================================
# PART B: CONDITIONAL LOGIT VS MULTINOMIAL LOGIT
# =============================================================================

# =============================================================================
# B: HELPER — Estimate multinomial logit on realised unions
# =============================================================================

estimate_mlogit_endogamy <- function(choice_data, agents) {
  #' Estimate multinomial logit on the realised unions.
  #'
  #' The outcome is the chosen partner's group. Predictors are the ego's
  #' group and education. Returns endogamy measures comparable to clogit.
  #'
  #' @param choice_data data.frame with chosen column
  #' @param agents data.frame with id, group, education
  #' @return list with mlogit coefficients and endogamy measures

  # Extract realised unions
  chosen <- choice_data[choice_data$chosen == 1, ]
  unions <- data.frame(
    ego_id = chosen$chooser_id,
    ego_group = chosen$ego_group,
    ego_edu = chosen$ego_education,
    partner_group = chosen$alt_group,
    partner_edu = chosen$alt_education
  )
  unions$partner_group_f <- factor(unions$partner_group)
  unions$ego_group_f <- factor(unions$ego_group)
  unions$same_group <- as.numeric(unions$ego_group == unions$partner_group)

  n_groups <- length(unique(agents$group))
  group_props <- as.numeric(table(agents$group) / nrow(agents))

  # --- Multinomial logit ---
  mfit <- tryCatch({
    multinom(partner_group_f ~ ego_group_f + ego_edu, data = unions, trace = FALSE)
  }, error = function(e) NULL)

  if (is.null(mfit)) {
    return(list(
      converged = FALSE,
      endogamy_predicted = rep(NA, n_groups),
      endogamy_expected = group_props,
      endogamy_ratio = rep(NA, n_groups),
      overall_endogamy_rate = mean(unions$same_group),
      same_group_coef = NA,
      mlogit_coefs = NA
    ))
  }

  # Predicted same-group probability for a representative ego from each group
  endogamy_predicted <- numeric(n_groups)
  for (g in 1:n_groups) {
    newdata <- data.frame(
      ego_group_f = factor(g, levels = levels(unions$ego_group_f)),
      ego_edu = mean(unions$ego_edu[unions$ego_group == g])
    )
    pred <- predict(mfit, newdata = newdata, type = "probs")
    if (is.matrix(pred)) pred <- pred[1, ]
    endogamy_predicted[g] <- pred[g]
  }

  # Endogamy ratio: predicted / expected under random matching
  endogamy_ratio <- endogamy_predicted / group_props

  # --- Simple logistic regression for direct endogamy coefficient ---
  # Controls for group sizes via offset
  unions$log_odds_expected <- log(group_props[unions$ego_group] /
                                    (1 - group_props[unions$ego_group]))
  gfit <- tryCatch({
    glm(same_group ~ ego_edu + ego_group_f + offset(log_odds_expected),
        family = binomial, data = unions)
  }, error = function(e) NULL)

  same_group_coef <- if (!is.null(gfit)) coef(gfit)["(Intercept)"] else NA

  list(
    converged = TRUE,
    endogamy_predicted = endogamy_predicted,
    endogamy_expected = group_props,
    endogamy_ratio = endogamy_ratio,
    mean_endogamy_ratio = mean(endogamy_ratio),
    overall_endogamy_rate = mean(unions$same_group),
    same_group_coef = same_group_coef,
    mlogit_coefs = if (!is.null(mfit)) coef(mfit) else NA,
    logistic_model = gfit
  )
}


# =============================================================================
# B: SIMULATION FUNCTION — Both clogit and mlogit on same data
# =============================================================================

run_comparison_sim <- function(
    N = N_AGENTS,
    J = J_ALTS,
    group_proportions = GROUP_PROPS,
    edu_means = EDU_MEANS,
    edu_sds = EDU_SDS,
    beta_edu = BETA_HOMOPHILY,
    beta_sg = BETA_SAME_GROUP,
    # Matching type: "one_sided" or "bilateral"
    matching_type = "one_sided",
    bilateral_threshold = -0.5,
    bilateral_rounds = 20,
    # Local market structure (Part B4-B5)
    use_local_markets = FALSE,
    n_markets = 5,
    market_segregation = 2.0,
    verbose = FALSE
) {
  n_groups <- length(group_proportions)

  # 1. Generate agents
  agents <- generate_agents(
    N = N, n_groups = n_groups,
    group_proportions = group_proportions,
    edu_means = edu_means, edu_sds = edu_sds
  )

  # 2. Optionally assign to local markets
  if (use_local_markets) {
    agents$market <- NA_integer_
    market_group_props <- matrix(NA, n_markets, n_groups)

    for (m in 1:n_markets) {
      focus <- ((m - 1) %% n_groups) + 1
      alpha <- rep(1, n_groups)
      alpha[focus] <- market_segregation
      market_group_props[m, ] <- alpha / sum(alpha)
    }

    for (i in 1:N) {
      g <- agents$group[i]
      market_weights <- market_group_props[, g]
      market_weights <- market_weights / sum(market_weights)
      agents$market[i] <- sample(1:n_markets, 1, prob = market_weights)
    }
  }

  # 3. Generate unions and construct choice sets
  if (matching_type == "bilateral") {
    # --- Two-sided matching ---
    bm <- bilateral_matching(
      agents, beta = c(beta_edu, beta_sg),
      var_names = c("edu_diff", "same_group"),
      n_rounds = bilateral_rounds, threshold = bilateral_threshold,
      verbose = FALSE
    )

    if (nrow(bm$matches) < 20) {
      warning("Bilateral matching produced too few pairs")
      return(NULL)
    }

    # Construct choice sets from matched pairs
    if (use_local_markets) {
      # For local markets + bilateral: choice sets drawn from local market
      all_ids <- agents$id
      cs_list <- vector("list", 2 * nrow(bm$matches))

      for (m_idx in 1:nrow(bm$matches)) {
        ego_a <- bm$matches$ego_id[m_idx]
        ego_b <- bm$matches$alt_id[m_idx]

        for (direction in 1:2) {
          chooser <- if (direction == 1) ego_a else ego_b
          partner <- if (direction == 1) ego_b else ego_a

          # Sample from local market
          market_i <- agents$market[agents$id == chooser]
          local_ids <- agents$id[agents$market == market_i &
                                   agents$id != chooser &
                                   agents$id != partner]
          j_actual <- min(J - 1, length(local_ids))
          if (j_actual < 2) next

          alts <- c(sample(local_ids, size = j_actual, replace = FALSE), partner)
          ego_row <- agents[agents$id == chooser, , drop = FALSE]
          alt_rows <- agents[match(alts, agents$id), , drop = FALSE]

          cs <- data.frame(chooser_id = chooser, alt_id = alts)
          for (col in setdiff(names(agents), "id")) {
            cs[[paste0("ego_", col)]] <- ego_row[[col]]
            cs[[paste0("alt_", col)]] <- alt_rows[[col]]
          }
          cs$chosen <- as.numeric(cs$alt_id == partner)

          list_idx <- (m_idx - 1) * 2 + direction
          cs_list[[list_idx]] <- cs
        }
      }
      cd <- bind_rows(cs_list[!sapply(cs_list, is.null)])
    } else {
      cd <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    }
    cd <- compute_match_variables(cd)

  } else {
    # --- One-sided matching (original) ---
    if (use_local_markets) {
      all_ids <- agents$id
      cs_list <- vector("list", N)

      for (idx in 1:N) {
        i <- agents$id[idx]
        market_i <- agents$market[idx]
        local_ids <- agents$id[agents$market == market_i & agents$id != i]

        n_available <- length(local_ids)
        j_actual <- min(J, n_available)
        if (j_actual < 2) next

        alts <- sample(local_ids, size = j_actual, replace = FALSE)
        ego <- agents[agents$id == i, , drop = FALSE]
        alt_df <- agents[agents$id %in% alts, , drop = FALSE]

        cs <- data.frame(chooser_id = i, alt_id = alt_df$id)
        for (col in setdiff(names(agents), "id")) {
          cs[[paste0("ego_", col)]] <- ego[[col]]
          cs[[paste0("alt_", col)]] <- alt_df[[col]]
        }
        cs_list[[idx]] <- cs
      }
      cd <- bind_rows(cs_list[!sapply(cs_list, is.null)])
    } else {
      cd <- construct_choice_sets(agents, J = J)
    }
    cd <- compute_match_variables(cd)
    cd <- simulate_choices(cd, beta = c(beta_edu, beta_sg),
                           var_names = c("edu_diff", "same_group"))
  }

  # 4. Diagnostics
  if (verbose) {
    cat(sprintf("\n  Matching type: %s\n", matching_type))
    describe_population(agents)
    if (use_local_markets) {
      cat("\n  ---- LOCAL MARKETS ----\n")
      for (m in 1:n_markets) {
        market_agents <- agents[agents$market == m, ]
        grp_counts <- table(market_agents$group)
        grp_pcts <- round(100 * grp_counts / nrow(market_agents))
        cat(sprintf("    Market %d: n = %d, groups: %s\n",
                    m, nrow(market_agents),
                    paste(sprintf("g%d=%s%%", as.integer(names(grp_pcts)), grp_pcts),
                          collapse = ", ")))
      }
    }
    describe_choice_sets(cd)
    describe_unions(cd)
  }

  # 5. Estimate clogit
  clogit_res <- estimate_clogit(cd, "edu_diff + same_group")

  # 6. Estimate multinomial logit on realised unions
  mlogit_res <- estimate_mlogit_endogamy(cd, agents)

  # 7. If local markets, also estimate mlogit with market controls
  mlogit_controlled <- NULL
  if (use_local_markets) {
    chosen <- cd[cd$chosen == 1, ]
    unions_ctrl <- data.frame(
      ego_group = chosen$ego_group,
      ego_edu = chosen$ego_education,
      partner_group = factor(chosen$alt_group),
      ego_group_f = factor(chosen$ego_group),
      ego_market = chosen$ego_market,
      same_group = as.numeric(chosen$ego_group == chosen$alt_group)
    )

    if ("ego_market" %in% names(unions_ctrl)) {
      for (g in 1:n_groups) {
        local_prop <- tapply(agents$group == g, agents$market, mean)
        unions_ctrl[[paste0("local_prop_g", g)]] <-
          local_prop[as.character(unions_ctrl$ego_market)]
      }

      ctrl_formula <- as.formula(paste(
        "partner_group ~ ego_group_f + ego_edu +",
        paste0("local_prop_g", 1:n_groups, collapse = " + ")
      ))
      mfit_ctrl <- tryCatch({
        multinom(ctrl_formula, data = unions_ctrl, trace = FALSE)
      }, error = function(e) NULL)

      if (!is.null(mfit_ctrl)) {
        endogamy_pred_ctrl <- numeric(n_groups)
        group_props <- as.numeric(table(agents$group) / nrow(agents))
        for (g in 1:n_groups) {
          newdata <- unions_ctrl[unions_ctrl$ego_group == g, ][1, , drop = FALSE]
          newdata$ego_group_f <- factor(g, levels = levels(unions_ctrl$ego_group_f))
          newdata$ego_edu <- mean(unions_ctrl$ego_edu[unions_ctrl$ego_group == g])
          for (gg in 1:n_groups) {
            newdata[[paste0("local_prop_g", gg)]] <- group_props[gg]
          }
          pred <- predict(mfit_ctrl, newdata = newdata, type = "probs")
          if (is.matrix(pred)) pred <- pred[1, ]
          endogamy_pred_ctrl[g] <- pred[g]
        }
        mlogit_controlled <- list(
          endogamy_predicted = endogamy_pred_ctrl,
          endogamy_ratio = endogamy_pred_ctrl / group_props,
          mean_endogamy_ratio = mean(endogamy_pred_ctrl / group_props)
        )
      }
    }
  }

  list(
    clogit = clogit_res,
    mlogit = mlogit_res,
    mlogit_controlled = mlogit_controlled
  )
}


# =============================================================================
# B: RUN COMPARISON TESTS
# =============================================================================

run_part_B <- function() {
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("PART B: Conditional Logit vs Multinomial Logit\n")
  cat("=", rep("=", 70), "\n")

  true_clogit <- c(edu_diff = BETA_HOMOPHILY, same_group = BETA_SAME_GROUP)
  # Use oracle (one-sided) only — bilateral matching tested in Module 4
  matching_types <- c("one_sided")
  all_B <- list()

  # Experiment definitions
  experiments_B <- list(
    B1 = list(
      name = "B1: Balanced groups (33/33/33)",
      group_proportions = c(1/3, 1/3, 1/3),
      edu_means = c(12, 12, 12), edu_sds = c(2, 2, 2),
      use_local_markets = FALSE
    ),
    B2 = list(
      name = "B2: Unbalanced (70/20/10), different edu",
      group_proportions = GROUP_PROPS,
      edu_means = EDU_MEANS, edu_sds = EDU_SDS,
      use_local_markets = FALSE
    ),
    B3 = list(
      name = "B3: Extreme (85/10/5), different edu",
      group_proportions = c(0.85, 0.10, 0.05),
      edu_means = c(12, 9, 15), edu_sds = c(2, 3, 1),
      use_local_markets = FALSE
    ),
    B4 = list(
      name = "B4: Local markets (segregation=3)",
      group_proportions = GROUP_PROPS,
      edu_means = EDU_MEANS, edu_sds = EDU_SDS,
      use_local_markets = TRUE, n_markets = 5, market_segregation = 3.0
    ),
    B5 = list(
      name = "B5: Local markets (segregation=10)",
      group_proportions = GROUP_PROPS,
      edu_means = EDU_MEANS, edu_sds = EDU_SDS,
      use_local_markets = TRUE, n_markets = 5, market_segregation = 10.0
    )
  )

  for (mtype in matching_types) {
    cat(sprintf("\n\n========== MATCHING TYPE: %s ==========\n",
                toupper(mtype)))

    for (exp_name in names(experiments_B)) {
      exp <- experiments_B[[exp_name]]
      label <- paste(exp_name, mtype, sep = "_")
      cat(sprintf("\n\n--- %s (%s) ---\n", exp$name, mtype))
      if (exp$use_local_markets) {
        cat("  Clogit conditions on local choice sets; mlogit does not.\n")
      }

      # Diagnostic run
      set.seed(99)
      extra_args <- list(
        group_proportions = exp$group_proportions,
        edu_means = exp$edu_means, edu_sds = exp$edu_sds,
        matching_type = mtype,
        use_local_markets = exp$use_local_markets,
        verbose = TRUE
      )
      if (exp$use_local_markets) {
        extra_args$n_markets <- exp$n_markets
        extra_args$market_segregation <- exp$market_segregation
      }
      do.call(run_comparison_sim, extra_args)

      # MC runs
      mc_args <- list(
        sim_fn = run_comparison_sim, R = R_SIMS, n_cores = N_CORES, seed = 42,
        group_proportions = exp$group_proportions,
        edu_means = exp$edu_means, edu_sds = exp$edu_sds,
        matching_type = mtype,
        use_local_markets = exp$use_local_markets
      )
      if (exp$use_local_markets) {
        mc_args$n_markets <- exp$n_markets
        mc_args$market_segregation <- exp$market_segregation
      }
      mc <- do.call(run_monte_carlo, mc_args)

      # Clogit recovery
      rec <- assess_recovery(mc, true_clogit, "clogit")
      cat(sprintf("\n  Clogit recovery (%s):\n", mtype))
      print(rec[, c("parameter", "true_value", "mean_estimate", "bias",
                     "rel_bias_pct", "coverage_95")])

      # Mlogit endogamy ratios
      mlogit_ratios <- sapply(mc, function(x) x$mlogit$mean_endogamy_ratio)
      mlogit_rates <- sapply(mc, function(x) x$mlogit$overall_endogamy_rate)
      cat(sprintf("\n  Mlogit mean endogamy ratio: %.2f (sd=%.2f)\n",
                  mean(mlogit_ratios, na.rm = TRUE), sd(mlogit_ratios, na.rm = TRUE)))

      # Per-group ratios for unbalanced conditions
      n_groups <- length(exp$group_proportions)
      if (n_groups <= 5) {
        for (g in 1:n_groups) {
          ratios_g <- sapply(mc, function(x) x$mlogit$endogamy_ratio[g])
          cat(sprintf("    Group %d endogamy ratio: %.2f (sd=%.2f)\n",
                      g, mean(ratios_g, na.rm = TRUE), sd(ratios_g, na.rm = TRUE)))
        }
      }

      # Mlogit with controls (local markets only)
      if (exp$use_local_markets) {
        ctrl_ratios <- sapply(mc, function(x) {
          if (!is.null(x$mlogit_controlled))
            x$mlogit_controlled$mean_endogamy_ratio else NA
        })
        cat(sprintf("  Mlogit (with local controls): %.2f (sd=%.2f)\n",
                    mean(ctrl_ratios, na.rm = TRUE), sd(ctrl_ratios, na.rm = TRUE)))
      }

      all_B[[label]] <- mc
    }
  }

  # ---- Summary comparison ----
  cat("\n\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("SUMMARY: One-sided vs Bilateral — Clogit bias across conditions\n")
  cat("=", rep("=", 70), "\n\n")

  for (exp_name in names(experiments_B)) {
    cat(sprintf("--- %s ---\n", experiments_B[[exp_name]]$name))
    for (mtype in matching_types) {
      label <- paste(exp_name, mtype, sep = "_")
      mc <- all_B[[label]]
      if (is.null(mc)) next

      rec <- assess_recovery(mc, true_clogit, "clogit")
      cat(sprintf("  %s: ", mtype))
      for (i in 1:nrow(rec)) {
        cat(sprintf("%s: bias=%.3f (%.1f%%)  ",
                    rec$parameter[i], rec$bias[i], rec$rel_bias_pct[i]))
      }
      cat("\n")
    }
    cat("\n")
  }

  return(all_B)
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("=", rep("=", 70), "\n", sep = "")
cat("MODULE 1b: Identification Tests\n")
cat("=", rep("=", 70), "\n")
cat(sprintf("Settings: N=%d, J=%d, R=%d reps per condition\n",
            N_AGENTS, J_ALTS, R_SIMS))

dir.create("results", showWarnings = FALSE)

part_A <- run_part_A()
part_B <- run_part_B()

saveRDS(list(part_A = part_A, part_B = part_B),
        "results/module1b_identification_results.rds")

cat("\n\nAll results saved to results/\n")
cat("--- Module 1b complete ---\n")
