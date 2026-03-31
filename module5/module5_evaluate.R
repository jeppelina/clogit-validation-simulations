# =============================================================================
# Module 5 — Evaluation Functions
# =============================================================================
#
# Recovery metrics, cluster comparison, and diagnostics.
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Recovery metrics
# -----------------------------------------------------------------------------

compute_recovery <- function(est_matrix, true_matrix, K, group_names = NULL) {
  #' Compute recovery metrics comparing estimated to true K×K matrix.
  #'
  #' @param est_matrix K×K matrix of estimated coefficients
  #' @param true_matrix K×K matrix of true alpha_gh values
  #' @param K Number of groups
  #' @param group_names Character vector of group names (optional)
  #' @return List of recovery metrics

  flat_est  <- as.vector(est_matrix)
  flat_true <- as.vector(true_matrix)

  # Symmetric version of true: (alpha_gh + alpha_hg) / 2
  sym_true <- (true_matrix + t(true_matrix)) / 2
  flat_sym <- as.vector(sym_true)

  # Min version: min(alpha_gh, alpha_hg) — the bottleneck
  min_true <- pmin(true_matrix, t(true_matrix))
  flat_min <- as.vector(min_true)

  # Full matrix recovery
  r_full     <- cor(flat_est, flat_true, use = "complete.obs")
  rho_full   <- cor(flat_est, flat_true, method = "spearman", use = "complete.obs")
  r_sym      <- cor(flat_est, flat_sym, use = "complete.obs")
  r_min      <- cor(flat_est, flat_min, use = "complete.obs")

  # Diagonal only (endogamy gradient)
  diag_est  <- diag(est_matrix)
  diag_true <- diag(true_matrix)
  r_diag    <- cor(diag_est, diag_true, use = "complete.obs")
  rho_diag  <- cor(diag_est, diag_true, method = "spearman", use = "complete.obs")

  # Off-diagonal only
  mask <- !diag(rep(TRUE, K))
  off_est  <- est_matrix[mask]
  off_true <- true_matrix[mask]
  off_sym  <- sym_true[mask]
  r_off     <- cor(off_est, off_true, use = "complete.obs")
  rho_off   <- cor(off_est, off_true, method = "spearman", use = "complete.obs")
  r_off_sym <- cor(off_est, off_sym, use = "complete.obs")

  # RMSE
  rmse_full <- sqrt(mean((flat_est - flat_true)^2, na.rm = TRUE))
  rmse_diag <- sqrt(mean((diag_est - diag_true)^2, na.rm = TRUE))

  list(
    r_full      = r_full,
    rho_full    = rho_full,
    r_sym       = r_sym,
    r_min       = r_min,
    r_diag      = r_diag,
    rho_diag    = rho_diag,
    r_off       = r_off,
    rho_off     = rho_off,
    r_off_sym   = r_off_sym,
    rmse_full   = rmse_full,
    rmse_diag   = rmse_diag
  )
}


# -----------------------------------------------------------------------------
# 2. Asymmetry diagnostic
# -----------------------------------------------------------------------------

asymmetry_diagnostic <- function(est_matrix, true_matrix, K, group_names) {
  #' For each off-diagonal pair, check whether the estimated asymmetry

  #' (est_gh - est_hg) has the same sign as the true asymmetry (alpha_gh - alpha_hg).
  #'
  #' @return data.frame with pair, true_gh, true_hg, est_gh, est_hg,
  #'         true_asym, est_asym, sign_match

  results <- list()
  for (g in 1:(K-1)) {
    for (h in (g+1):K) {
      true_asym <- true_matrix[g, h] - true_matrix[h, g]
      est_asym  <- est_matrix[g, h] - est_matrix[h, g]
      results[[length(results) + 1]] <- data.frame(
        pair       = paste0(group_names[g], "-", group_names[h]),
        true_gh    = true_matrix[g, h],
        true_hg    = true_matrix[h, g],
        est_gh     = est_matrix[g, h],
        est_hg     = est_matrix[h, g],
        true_asym  = true_asym,
        est_asym   = est_asym,
        sign_match = sign(true_asym) == sign(est_asym),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results)
}


# -----------------------------------------------------------------------------
# 3. Cluster recovery (requires igraph)
# -----------------------------------------------------------------------------

recover_clusters <- function(coef_matrix, K, group_names = NULL) {
  #' Apply Walktrap community detection to a K×K coefficient matrix.
  #'
  #' Only positive coefficients are used as edge weights (as in Paper 2).
  #'
  #' @return List with $membership (integer vector), $modularity, $n_clusters

  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("igraph not installed. Cluster recovery skipped.")
    return(list(membership = rep(NA, K), modularity = NA, n_clusters = NA))
  }

  # Use only positive coefficients as edges
  adj <- coef_matrix
  adj[adj < 0] <- 0
  diag(adj) <- 0  # remove self-loops

  # Create directed weighted graph
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE)

  # Walktrap
  wt <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight, steps = 4)
  mem <- igraph::membership(wt)
  mod <- igraph::modularity(wt)

  if (!is.null(group_names)) names(mem) <- group_names

  list(membership = mem, modularity = mod, n_clusters = length(unique(mem)))
}


compare_clusters <- function(mem_true, mem_est) {
  #' Compare two cluster assignments using Adjusted Rand Index.
  #'
  #' If aricode is available, use it. Otherwise, compute manually.
  #' @return Scalar ARI value (1 = perfect, 0 = random, <0 = worse than random)

  if (requireNamespace("aricode", quietly = TRUE)) {
    return(aricode::ARI(mem_true, mem_est))
  }

  # Manual ARI computation
  n <- length(mem_true)
  tab <- table(mem_true, mem_est)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2))
  c <- sum(choose(colSums(tab), 2))
  d <- choose(n, 2)
  expected <- b * c / d
  max_val <- (b + c) / 2
  if (max_val == expected) return(1)
  (a - expected) / (max_val - expected)
}


# -----------------------------------------------------------------------------
# 4. Summary printer
# -----------------------------------------------------------------------------

print_recovery_comparison <- function(results_list, alpha, group_names) {
  #' Print a formatted comparison of recovery metrics across methods.
  #'
  #' @param results_list Named list of recovery metric lists
  #' @param alpha True alpha matrix

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("RECOVERY COMPARISON\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Header
  cat(sprintf("%-25s %8s %8s %8s %8s %8s %8s\n",
              "Method", "r_full", "ρ_full", "r_diag", "ρ_diag", "r_off", "ρ_off"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (name in names(results_list)) {
    r <- results_list[[name]]
    cat(sprintf("%-25s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                name,
                r$r_full, r$rho_full,
                r$r_diag, r$rho_diag,
                r$r_off, r$rho_off))
  }

  cat("\n")
  cat(sprintf("%-25s %8s %8s %8s\n", "Method", "vs_sym", "vs_min", "RMSE"))
  cat(paste(rep("-", 55), collapse = ""), "\n")
  for (name in names(results_list)) {
    r <- results_list[[name]]
    cat(sprintf("%-25s %8.3f %8.3f %8.3f\n",
                name, r$r_sym, r$r_min, r$rmse_full))
  }
}


cat("Module 5 evaluation functions loaded.\n")
