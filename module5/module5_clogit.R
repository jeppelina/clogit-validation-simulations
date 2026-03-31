# =============================================================================
# Module 5 — Clogit Pipeline (v2 — separate models per ego group)
# =============================================================================
#
# Estimates K separate clogit models, one per ego group.
# Each model has K partner-group dummies + log-distance.
# This mirrors Paper 1's approach and avoids overparameterization.
#
# =============================================================================

library(survival)

# -----------------------------------------------------------------------------
# 1. Construct choice sets
# -----------------------------------------------------------------------------

construct_choice_sets_m5 <- function(unions, men, women, J = 50) {
  #' For each observed union (men's perspective), build a choice set:
  #' actual partner + (J-1) random alternative women.

  N_unions <- nrow(unions)
  N_w <- nrow(women)

  rows <- vector("list", N_unions)

  for (u in 1:N_unions) {
    mi <- unions$man_id[u]
    wi <- unions$woman_id[u]
    mg <- unions$man_group[u]

    man_row <- which(men$id == mi)[1]
    man_x <- men$x[man_row]

    partner_row <- which(women$id == wi)[1]

    # Sample J-1 alternatives (exclude actual partner)
    alt_pool <- setdiff(seq_len(N_w), partner_row)
    n_alts <- min(J - 1, length(alt_pool))
    alt_rows <- sample(alt_pool, n_alts, replace = FALSE)

    all_rows <- c(partner_row, alt_rows)
    n_total <- length(all_rows)

    rows[[u]] <- data.frame(
      ego_id    = rep(mi, n_total),
      ego_group = rep(mg, n_total),
      alt_id    = women$id[all_rows],
      alt_group = women$group[all_rows],
      log_dist  = log(abs(man_x - women$x[all_rows]) + 0.5),
      chosen    = c(1L, rep(0L, n_total - 1)),
      strata_id = rep(u, n_total),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# 2. Estimate K separate clogit models
# -----------------------------------------------------------------------------

estimate_clogit_m5 <- function(choice_data, K, group_names = NULL) {
  #' Run one clogit per ego group. Each model estimates K partner-group
  #' dummies (one per alt_group) plus log_dist.
  #'
  #' Returns a K×K coefficient matrix where row g = ego group,
  #' column h = coefficient for partner group h from model g.
  #' The reference within each model is normalized so that the
  #' mean coefficient is 0 (centering).

  coef_matrix <- matrix(NA_real_, nrow = K, ncol = K)
  delta_vec <- rep(NA_real_, K)
  converged_vec <- rep(FALSE, K)

  for (g in 1:K) {
    # Subset to egos from group g
    sub <- choice_data[choice_data$ego_group == g, ]

    if (nrow(sub) == 0 || sum(sub$chosen) < 3) {
      # Too few unions from this group
      next
    }

    # Create partner-group dummies
    for (h in 1:K) {
      sub[[paste0("pg_", h)]] <- as.integer(sub$alt_group == h)
    }

    # Build formula: use all K dummies (one will be dropped by clogit for
    # identification, or we can drop the largest group for stability)
    # Actually, include all K — clogit will handle rank deficiency.
    # Safer: drop one explicitly. Drop the group with most alternatives.
    alt_group_counts <- table(sub$alt_group)
    ref_h <- as.integer(names(which.max(alt_group_counts)))

    pg_vars <- paste0("pg_", setdiff(1:K, ref_h))
    fml_rhs <- paste(c(pg_vars, "log_dist"), collapse = " + ")

    # clogit accepts chosen ~ ... + strata(id) directly (no Surv wrapper needed)
    fml <- as.formula(paste("chosen ~", fml_rhs, "+ strata(strata_id)"))

    fit <- tryCatch(
      clogit(fml, data = sub, method = "efron"),
      error = function(e) {
        if (!is.null(group_names)) {
          cat(sprintf("    clogit failed for group %s: %s\n", group_names[g], e$message))
        }
        NULL
      },
      warning = function(w) {
        # Catch convergence warnings but still try to get coefficients
        suppressWarnings(clogit(fml, data = sub, method = "efron"))
      }
    )

    if (is.null(fit)) next

    converged_vec[g] <- TRUE
    coefs <- coef(fit)

    # Fill in the coefficient row
    for (h in 1:K) {
      varname <- paste0("pg_", h)
      if (varname %in% names(coefs)) {
        coef_matrix[g, h] <- coefs[varname]
      } else if (h == ref_h) {
        coef_matrix[g, h] <- 0  # reference category
      }
    }

    if ("log_dist" %in% names(coefs)) {
      delta_vec[g] <- coefs["log_dist"]
    }
  }

  # Center each row so mean = 0 (makes coefficients comparable across groups)
  for (g in 1:K) {
    if (!any(is.na(coef_matrix[g, ]))) {
      coef_matrix[g, ] <- coef_matrix[g, ] - mean(coef_matrix[g, ])
    }
  }

  n_converged <- sum(converged_vec)

  list(
    coef_matrix = coef_matrix,
    delta_vec   = delta_vec,
    converged   = converged_vec,
    n_converged = n_converged
  )
}


cat("Module 5 clogit pipeline loaded (v2: separate models per ego group).\n")
