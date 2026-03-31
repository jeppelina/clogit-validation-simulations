# =============================================================================
# Module 5 — Data-Generating Process
# =============================================================================
#
# Functions for:
#   1. Creating populations with K groups and 1D spatial positions
#   2. Computing pairwise utility from the K×K alpha matrix
#   3. Bilateral matching (gravity-based, with mutual acceptance)
#   4. One-sided matching (gravity-based, for sanity check only)
#
# =============================================================================

library(survival)

# -----------------------------------------------------------------------------
# 1. Create population
# -----------------------------------------------------------------------------

create_population <- function(group_sizes, geo_centers, geo_spread, group_names) {
  #' Generate men and women with group membership and spatial positions.
  #'
  #' @param group_sizes Integer vector of size K (per gender)
  #' @param geo_centers Named numeric vector of geographic centers
  #' @param geo_spread Scalar: within-group SD of positions
  #' @param group_names Character vector of group names
  #' @return List with $men and $women data.frames (id, group, group_name, x)

  K <- length(group_sizes)
  N_total <- sum(group_sizes)

  make_agents <- function(prefix) {
    dfs <- list()
    cum_id <- 0
    for (g in 1:K) {
      n_g <- group_sizes[g]
      center <- geo_centers[g]
      dfs[[g]] <- data.frame(
        id     = (cum_id + 1):(cum_id + n_g),
        group  = g,
        group_name = group_names[g],
        x      = rnorm(n_g, mean = center, sd = geo_spread),
        stringsAsFactors = FALSE
      )
      cum_id <- cum_id + n_g
    }
    do.call(rbind, dfs)
  }

  list(men = make_agents("m"), women = make_agents("w"))
}


# -----------------------------------------------------------------------------
# 2. Pairwise utility
# -----------------------------------------------------------------------------

compute_utility <- function(man, woman, alpha, delta) {
  #' Deterministic utility of (man, woman) pair.
  #'
  #' V = alpha[man_group, woman_group] - delta * |x_man - x_woman|
  #'
  #' @param man  One-row data.frame (or list) with $group and $x

  #' @param woman One-row data.frame (or list) with $group and $x
  #' @param alpha K×K preference matrix
  #' @param delta Distance weight (scalar, positive = distance costly)
  #' @return Scalar deterministic utility V

  alpha[man$group, woman$group] - delta * abs(man$x - woman$x)
}

# Vectorized version for speed
compute_utility_vec <- function(man_groups, man_x, woman_groups, woman_x,
                                alpha, delta) {
  #' Vectorized utility for vectors of pairs.
  #' @return Numeric vector of deterministic utilities
  alpha[cbind(man_groups, woman_groups)] - delta * abs(man_x - woman_x)
}


# -----------------------------------------------------------------------------
# 3. Bilateral matching (gravity-based with mutual threshold)
# -----------------------------------------------------------------------------

bilateral_gravity_match <- function(men, women, alpha, delta,
                                    acceptance_rate = 0.50,
                                    n_rounds = 80,
                                    n_calibration = 5000,
                                    verbose = FALSE) {
  #' Bilateral matching: random encounters, both must have U > threshold.
  #'
  #' The threshold is calibrated so that the one-sided acceptance rate
  #' (P(U_one_side > threshold)) equals the target. Joint acceptance rate
  #' is approximately acceptance_rate^2 for independent pairs.
  #'
  #' Man's utility uses alpha[man_group, woman_group].
  #' Woman's utility uses alpha[woman_group, man_group] (transpose).
  #'
  #' @param men data.frame with id, group, x
  #' @param women data.frame with id, group, x
  #' @param alpha K×K preference matrix (row = ego's group)
  #' @param delta Distance weight
  #' @param acceptance_rate Target one-sided acceptance rate for calibration
  #' @param n_rounds Number of random encounter rounds
  #' @param n_calibration Number of sample pairs for threshold calibration
  #' @param verbose Print progress
  #' @return List with $unions (data.frame: man_id, woman_id, man_group, woman_group),
  #'         $pairing_rate, $threshold

  N_m <- nrow(men)
  N_w <- nrow(women)

  # --- Calibrate threshold ---
  samp_m <- sample(N_m, n_calibration, replace = TRUE)
  samp_w <- sample(N_w, n_calibration, replace = TRUE)
  V_samp <- compute_utility_vec(men$group[samp_m], men$x[samp_m],
                                women$group[samp_w], women$x[samp_w],
                                alpha, delta)
  eps_samp <- -log(-log(runif(n_calibration)))
  U_samp <- V_samp + eps_samp
  threshold <- as.numeric(quantile(U_samp, 1 - acceptance_rate))

  if (verbose) cat(sprintf("  Calibrated threshold: %.3f (target rate: %.0f%%)\n",
                           threshold, acceptance_rate * 100))

  # --- Run matching ---
  matched_m <- logical(N_m)
  matched_w <- logical(N_w)
  unions <- list()
  n_unions <- 0

  for (round in 1:n_rounds) {
    avail_m <- which(!matched_m)
    avail_w <- which(!matched_w)
    if (length(avail_m) < 2 || length(avail_w) < 2) break

    # Random encounters: pair each available man with a random available woman
    n_enc <- min(length(avail_m), length(avail_w))
    enc_m <- sample(avail_m, n_enc)
    enc_w <- sample(avail_w, n_enc)

    # Compute utilities for all encounters
    V_mw <- compute_utility_vec(men$group[enc_m], men$x[enc_m],
                                women$group[enc_w], women$x[enc_w],
                                alpha, delta)
    V_wm <- compute_utility_vec(women$group[enc_w], women$x[enc_w],
                                men$group[enc_m], men$x[enc_m],
                                alpha, delta)

    eps_m <- -log(-log(runif(n_enc)))
    eps_w <- -log(-log(runif(n_enc)))

    U_m <- V_mw + eps_m
    U_w <- V_wm + eps_w

    # Bilateral acceptance
    accept <- (U_m > threshold) & (U_w > threshold)

    if (any(accept)) {
      acc_idx <- which(accept)
      for (idx in acc_idx) {
        mi <- enc_m[idx]
        wi <- enc_w[idx]
        if (!matched_m[mi] && !matched_w[wi]) {
          matched_m[mi] <- TRUE
          matched_w[wi] <- TRUE
          n_unions <- n_unions + 1
          unions[[n_unions]] <- c(men$id[mi], women$id[wi],
                                  men$group[mi], women$group[wi])
        }
      }
    }

    if (verbose && round %% 20 == 0) {
      cat(sprintf("    Round %d: %d unions, %.0f%% paired\n",
                  round, n_unions, 100 * 2 * n_unions / (N_m + N_w)))
    }
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

  pairing_rate <- 2 * n_unions / (N_m + N_w)

  if (verbose) cat(sprintf("  Bilateral done: %d unions (%.0f%% paired)\n",
                           n_unions, 100 * pairing_rate))

  list(unions = union_df, pairing_rate = pairing_rate, threshold = threshold)
}


# -----------------------------------------------------------------------------
# 4. One-sided gravity matching (sanity check only)
# -----------------------------------------------------------------------------

onesided_gravity_match <- function(men, women, alpha, delta, verbose = FALSE) {
  #' One-sided sequential gravity matching with pool depletion.
  #'
  #' Each man probabilistically selects an available woman.
  #' P(choose j) ∝ exp(alpha[g_i, g_j] - delta * |x_i - x_j|)
  #' No mutual acceptance — chosen partner always accepts.
  #'
  #' @return List with $unions data.frame, $pairing_rate

  N_m <- nrow(men)
  N_w <- nrow(women)
  available <- rep(TRUE, N_w)
  unions <- list()
  n_unions <- 0

  order <- sample(N_m)

  for (i in order) {
    avail_idx <- which(available)
    if (length(avail_idx) == 0) break

    # Compute weights for all available women
    V <- compute_utility_vec(
      rep(men$group[i], length(avail_idx)), rep(men$x[i], length(avail_idx)),
      women$group[avail_idx], women$x[avail_idx],
      alpha, delta
    )
    # Add Gumbel noise
    eps <- -log(-log(runif(length(avail_idx))))
    U <- V + eps

    # Pick the highest
    best <- which.max(U)
    j <- avail_idx[best]

    available[j] <- FALSE
    n_unions <- n_unions + 1
    unions[[n_unions]] <- c(men$id[i], women$id[j], men$group[i], women$group[j])
  }

  union_mat <- do.call(rbind, unions)
  union_df <- data.frame(
    man_id     = union_mat[, 1],
    woman_id   = union_mat[, 2],
    man_group  = union_mat[, 3],
    woman_group = union_mat[, 4]
  )

  list(unions = union_df, pairing_rate = 2 * n_unions / (N_m + N_w))
}


# -----------------------------------------------------------------------------
# 5. Union matrix
# -----------------------------------------------------------------------------

unions_to_matrix <- function(unions, K) {
  #' Convert unions data.frame to K×K frequency matrix.
  #' Rows = man's group, Columns = woman's group.
  mat <- matrix(0, nrow = K, ncol = K)
  for (r in 1:nrow(unions)) {
    mat[unions$man_group[r], unions$woman_group[r]] <-
      mat[unions$man_group[r], unions$woman_group[r]] + 1
  }
  mat
}

cat("Module 5 DGP functions loaded.\n")
