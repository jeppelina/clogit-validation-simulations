# =============================================================================
# Module 5 — Paper 2 Pipeline
# =============================================================================
#
# Implements Paper 2's method:
#   1. Null simulation (geography-only rematching with pool depletion)
#   2. Stacking observed + null into 3D table
#   3. Fitting saturated log-linear model
#   4. Extracting three-way interaction coefficients
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Null simulation (geography-only matching)
# -----------------------------------------------------------------------------

null_simulation <- function(men, women, delta,
                            search_radius_start = 4.0,
                            search_radius_max = 20.0,
                            search_radius_step = 2.0) {
  #' Paper 2's null model: rematch within geographic proximity, ignoring ancestry.
  #'
  #' For each man, find available women within expanding radius of his position.
  #' Select randomly (weighted by distance decay). No ancestry information used.
  #'
  #' @param men data.frame with id, group, x
  #' @param women data.frame with id, group, x
  #' @param delta Distance decay parameter
  #' @param search_radius_start Initial search radius
  #' @param search_radius_max Maximum search radius
  #' @param search_radius_step Step for expanding radius
  #' @return List with $unions data.frame, $n_matched

  N_m <- nrow(men)
  N_w <- nrow(women)
  available <- rep(TRUE, N_w)
  unions <- list()
  n_unions <- 0

  # Pre-compute women's x for speed
  w_x <- women$x
  w_group <- women$group

  order <- sample(N_m)

  for (i in order) {
    avail_idx <- which(available)
    if (length(avail_idx) == 0) break

    man_x <- men$x[i]

    # Find candidates within expanding radius
    candidates <- NULL
    radius <- search_radius_start

    while (is.null(candidates) || length(candidates) < 10) {
      dists <- abs(man_x - w_x[avail_idx])
      candidates <- avail_idx[dists <= radius]
      if (length(candidates) >= 10 || radius >= search_radius_max) break
      radius <- radius + search_radius_step
    }

    if (length(candidates) == 0) {
      # Fallback: all available
      candidates <- avail_idx
    }

    # Select with distance-decay weighting (no ancestry)
    if (length(candidates) == 1) {
      j <- candidates[1]
    } else {
      d_cand <- abs(man_x - w_x[candidates])
      weights <- exp(-delta * d_cand)
      weights <- weights / sum(weights)
      j <- candidates[sample.int(length(candidates), 1, prob = weights)]
    }

    available[j] <- FALSE
    n_unions <- n_unions + 1
    unions[[n_unions]] <- c(men$id[i], women$id[j],
                            men$group[i], w_group[j])
  }

  if (n_unions == 0) {
    union_df <- data.frame(man_id = integer(0), woman_id = integer(0),
                           man_group = integer(0), woman_group = integer(0))
  } else {
    union_mat <- do.call(rbind, unions)
    union_df <- data.frame(
      man_id     = union_mat[, 1],
      woman_id   = union_mat[, 2],
      man_group  = union_mat[, 3],
      woman_group = union_mat[, 4]
    )
  }

  list(unions = union_df, n_matched = n_unions)
}


# -----------------------------------------------------------------------------
# 2. Run multiple null simulations and average
# -----------------------------------------------------------------------------

run_null_sims <- function(men, women, delta, n_sims, K,
                          verbose = FALSE) {
  #' Run n_sims null simulations and return the average K×K matrix.
  #'
  #' @return List with $avg_null (K×K matrix), $null_mats (list of matrices),
  #'         $sd_null (K×K matrix of across-sim SDs)

  null_mats <- vector("list", n_sims)

  for (s in 1:n_sims) {
    sim <- null_simulation(men, women, delta)
    null_mats[[s]] <- unions_to_matrix(sim$unions, K)
    if (verbose && s %% 10 == 0) {
      cat(sprintf("    Null sim %d/%d\n", s, n_sims))
    }
  }

  # Average and SD across simulations
  mat_array <- array(unlist(null_mats), dim = c(K, K, n_sims))
  avg_null <- apply(mat_array, c(1, 2), mean)
  sd_null  <- apply(mat_array, c(1, 2), sd)

  list(avg_null = avg_null, null_mats = null_mats, sd_null = sd_null)
}


# -----------------------------------------------------------------------------
# 3. Log-linear model: three-way interaction
# -----------------------------------------------------------------------------

fit_loglinear <- function(obs_mat, null_mats, K, z = 0.5) {
  #' Fit a saturated log-linear model to the stacked observed + null tables.
  #'
  #' Replicates Paper 2's actual method: stack the observed K×K matrix with
  #' ALL S null simulation matrices into a K × K × (S+1) table, then fit
  #' the saturated log-linear model and extract the three-way interaction
  #' for the observed scenario (k=1).
  #'
  #' @param obs_mat K×K matrix of observed union counts
  #' @param null_mats LIST of K×K matrices from null simulations (length S)
  #' @param K Number of groups
  #' @param z Small constant added to zero cells
  #' @return List with $three_way (K×K matrix of three-way interaction for k=1),
  #'         $crude_ratio (K×K matrix of log(obs/mean_null)),
  #'         $loglinear_converged (logical)

  S <- length(null_mats)
  n_scenarios <- S + 1  # 1 observed + S nulls

  # --- Construct the 3D table: K × K × (S+1) ---
  tab <- array(0, dim = c(K, K, n_scenarios))
  tab[, , 1] <- obs_mat
  for (s in 1:S) {
    tab[, , s + 1] <- null_mats[[s]]
  }

  # Add z to all zero cells (not to all cells — only zeros)
  tab[tab == 0] <- z

  # --- Fit saturated log-linear model ---
  # Margins for saturated model: include the three-way interaction
  # Dim 1 = man's group, Dim 2 = woman's group, Dim 3 = scenario
  fit <- tryCatch({
    loglin(tab, margin = list(c(1, 2, 3)),
           fit = TRUE, param = TRUE, print = FALSE)
  }, error = function(e) {
    warning("loglin() failed: ", e$message)
    return(NULL)
  })

  converged <- !is.null(fit)

  # --- Extract three-way interaction for scenario k=1 (observed) ---
  if (converged && !is.null(fit$param)) {
    # fit$param contains all lambda parameters
    # The three-way interaction is in fit$param[[4]] (the 1.2.3 interaction)
    # It's a K × K × (S+1) array
    # We want the slice for k=1 (the observed scenario)
    three_way_full <- fit$param[["1.2.3"]]  # may also be indexed differently

    if (is.null(three_way_full)) {
      # Try numeric index — loglin names vary
      # The three-way is always the last element for a saturated model
      param_names <- names(fit$param)
      three_way_full <- fit$param[[length(fit$param)]]
    }

    if (!is.null(three_way_full) && length(dim(three_way_full)) == 3) {
      three_way <- three_way_full[, , 1]  # slice for observed scenario
    } else {
      warning("Could not extract three-way interaction from loglin params")
      three_way <- matrix(NA, K, K)
    }
  } else {
    # Fallback: manual extraction from the log-frequency table
    # The three-way interaction = log(F) with all lower-order terms removed
    # For a saturated model, F = observed counts
    log_F <- log(tab)

    # Remove all lower-order effects by iterated centering
    # Step 1: remove the three-way interaction by centering on each dimension
    # Actually, for extraction we want the OPPOSITE: remove lower-order, keep three-way
    #
    # Three-way interaction = log(F_ijk) centered across all three dimensions
    # = log(F_ijk) - mean_k(log(F_ijk)) - mean_j(log(F_ijk)) - mean_i(log(F_ijk))
    #   + mean_jk(log(F_ijk)) + mean_ik(log(F_ijk)) + mean_ij(log(F_ijk))
    #   - mean_ijk(log(F_ijk))

    m1 <- apply(log_F, 1, mean)  # mean over j,k for each i
    m2 <- apply(log_F, 2, mean)  # mean over i,k for each j
    m3 <- apply(log_F, 3, mean)  # mean over i,j for each k
    m12 <- apply(log_F, c(1,2), mean)  # mean over k for each i,j
    m13 <- apply(log_F, c(1,3), mean)  # mean over j for each i,k
    m23 <- apply(log_F, c(2,3), mean)  # mean over i for each j,k
    m123 <- mean(log_F)

    three_way_full <- array(0, dim = c(K, K, n_scenarios))
    for (i in 1:K) {
      for (j in 1:K) {
        for (k in 1:n_scenarios) {
          three_way_full[i,j,k] <- log_F[i,j,k] -
            m12[i,j] - m13[i,k] - m23[j,k] +
            m1[i] + m2[j] + m3[k] -
            m123
        }
      }
    }
    three_way <- three_way_full[, , 1]  # observed scenario
  }

  # --- Crude ratio (for comparison) ---
  avg_null <- Reduce("+", null_mats) / S
  crude_ratio <- log((obs_mat + z) / (avg_null + z))

  list(
    three_way          = three_way,
    crude_ratio        = crude_ratio,
    obs_mat            = obs_mat,
    avg_null           = avg_null,
    n_scenarios        = n_scenarios,
    loglinear_converged = converged,
    z                  = z
  )
}


cat("Module 5 Paper 2 pipeline loaded.\n")
