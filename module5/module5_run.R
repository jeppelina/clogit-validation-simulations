# =============================================================================
# Module 5 — Main Runner
# =============================================================================
#
# Monte Carlo experiment comparing Paper 2's simulation+loglinear pipeline
# against clogit with K² interactions for recovering known boundary structures.
#
# Core experiment: bilateral matching at 3 selectivity levels.
# Sanity check: one-sided matching (single run).
#
# Usage:
#   Rscript module5_run.R              # full run
#   Rscript module5_run.R quick        # quick test (5 MC reps)
#
# Output: results/module5_results.csv, results/module5_summary.md
#
# =============================================================================

cat("=" , rep("=", 69), "\n", sep = "")
cat("Module 5: Boundary Structure Recovery Validation\n")
cat(rep("=", 70), "\n\n", sep = "")

# --- Determine script directory robustly ---
# Works with Rscript, source(), and RStudio
get_script_dir <- function() {
  # Try commandArgs (Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  # Try sys.frame (source())
  for (i in sys.nframe():1) {
    if (!is.null(sys.frame(i)$ofile)) {
      return(dirname(normalizePath(sys.frame(i)$ofile)))
    }
  }
  # Fallback: working directory
  return(getwd())
}

script_dir <- get_script_dir()
cat(sprintf("Script directory: %s\n", script_dir))

# --- Load modules ---
source(file.path(script_dir, "module5_config.R"))
source(file.path(script_dir, "module5_dgp.R"))
source(file.path(script_dir, "module5_paper2.R"))
source(file.path(script_dir, "module5_clogit.R"))
source(file.path(script_dir, "module5_evaluate.R"))

# --- Quick mode? ---
args <- commandArgs(trailingOnly = TRUE)
if ("quick" %in% args) {
  N_MC <- 5
  N_NULL_SIMS <- 10
  J_CHOICE <- 30
  cat("*** QUICK MODE: 5 MC reps, 10 null sims ***\n\n")
}

# --- Results storage ---
results_file <- file.path(RESULTS_DIR, "module5_results.csv")

# Initialize results CSV
if (!file.exists(results_file)) {
  header <- data.frame(
    rep = integer(), condition = character(), acceptance_rate = numeric(),
    method = character(),
    r_full = numeric(), rho_full = numeric(),
    r_diag = numeric(), rho_diag = numeric(),
    r_off = numeric(), rho_off = numeric(),
    r_sym = numeric(), r_min = numeric(),
    rmse_full = numeric(), rmse_diag = numeric(),
    pairing_rate = numeric(),
    elapsed_sec = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(header, results_file, row.names = FALSE)
}

# --- Derived targets ---
ALPHA_SYM <- (ALPHA + t(ALPHA)) / 2
TRUE_CLUSTERS <- c(1, 1, 2, 3, 3)  # expected: {Nordic,European}, {LatAm}, {MENA,HornAfr}
names(TRUE_CLUSTERS) <- GROUP_NAMES

cat("True alpha matrix:\n")
print(round(ALPHA, 2))
cat("\nExpected clusters:", paste(GROUP_NAMES, "→", TRUE_CLUSTERS), "\n\n")


# =============================================================================
# MAIN LOOP
# =============================================================================

set.seed(SEED)

for (rep in 1:N_MC) {
  cat(sprintf("\n--- MC rep %d/%d ---\n", rep, N_MC))
  t0_rep <- Sys.time()

  # 1. Generate population
  pop <- create_population(GROUP_SIZES, GEO_CENTERS, GEO_SPREAD, GROUP_NAMES)
  men <- pop$men
  women <- pop$women

  # 2. Generate null simulations (shared across conditions for this rep)
  cat("  Running null simulations...\n")
  null_result <- run_null_sims(men, women, DELTA, N_NULL_SIMS, K, verbose = FALSE)
  avg_null <- null_result$avg_null

  # =========================================================================
  # BILATERAL CONDITIONS (the core experiment)
  # =========================================================================
  for (acc_rate in ACCEPTANCE_RATES) {
    cat(sprintf("  Bilateral (acceptance = %.0f%%)...\n", acc_rate * 100))
    t0 <- Sys.time()

    # Generate observed unions
    obs_result <- bilateral_gravity_match(men, women, ALPHA, DELTA,
                                           acceptance_rate = acc_rate,
                                           n_rounds = N_ROUNDS, verbose = FALSE)
    obs_unions <- obs_result$unions
    obs_mat <- unions_to_matrix(obs_unions, K)

    if (nrow(obs_unions) < 20) {
      cat("    WARNING: Very few unions formed. Skipping.\n")
      next
    }

    # --- Method 1: Paper 2 crude ratio ---
    p2_result <- fit_loglinear(obs_mat, null_result$null_mats, K)
    rec_crude <- compute_recovery(p2_result$crude_ratio, ALPHA, K, GROUP_NAMES)
    rec_3way  <- compute_recovery(p2_result$three_way, ALPHA, K, GROUP_NAMES)

    # --- Method 2: Clogit with K separate models ---
    cat("    Constructing choice sets...\n")
    cs <- construct_choice_sets_m5(obs_unions, men, women, J = J_CHOICE)
    cat("    Estimating clogit (K separate models)...\n")
    clogit_result <- estimate_clogit_m5(cs, K, GROUP_NAMES)
    cat(sprintf("    Clogit: %d/%d group models converged\n",
                clogit_result$n_converged, K))
    if (clogit_result$n_converged >= K - 1) {
      rec_clogit <- compute_recovery(clogit_result$coef_matrix, ALPHA, K, GROUP_NAMES)
    } else {
      rec_clogit <- list(r_full = NA, rho_full = NA, r_diag = NA, rho_diag = NA,
                          r_off = NA, rho_off = NA, r_sym = NA, r_min = NA,
                          rmse_full = NA, rmse_diag = NA)
    }

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # Save results
    for (method_info in list(
      list(name = "crude_ratio",    rec = rec_crude),
      list(name = "loglinear_3way", rec = rec_3way),
      list(name = "clogit_K2",     rec = rec_clogit)
    )) {
      row <- data.frame(
        rep = rep, condition = "bilateral", acceptance_rate = acc_rate,
        method = method_info$name,
        r_full = method_info$rec$r_full, rho_full = method_info$rec$rho_full,
        r_diag = method_info$rec$r_diag, rho_diag = method_info$rec$rho_diag,
        r_off = method_info$rec$r_off, rho_off = method_info$rec$rho_off,
        r_sym = method_info$rec$r_sym, r_min = method_info$rec$r_min,
        rmse_full = method_info$rec$rmse_full, rmse_diag = method_info$rec$rmse_diag,
        pairing_rate = obs_result$pairing_rate,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
      write.table(row, results_file, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    }

    cat(sprintf("    Done (%.0f sec). Pairing rate: %.0f%%\n",
                elapsed, 100 * obs_result$pairing_rate))
  }

  # =========================================================================
  # VALIDATION DIAGNOSTICS (first rep only)
  # =========================================================================
  if (rep == 1) {
    # --- Save validation diagnostics to a log file ---
    val_log <- file.path(RESULTS_DIR, "module5_validation.log")
    val_con <- file(val_log, open = "wt")
    # Helper: write to both console and log
    vcat <- function(...) { msg <- paste0(...); cat(msg); cat(msg, file = val_con) }
    vprint <- function(x) { capture.output(print(x)) |> paste(collapse = "\n") |>
                              (\(s) { cat(s, "\n"); cat(s, "\n", file = val_con) })() }

    vcat("\n  ======== VALIDATION DIAGNOSTICS (Rep 1) ========\n")

    # Pick the 50% bilateral condition for diagnostics
    diag_obs <- bilateral_gravity_match(men, women, ALPHA, DELTA,
                                         acceptance_rate = 0.50,
                                         n_rounds = N_ROUNDS, verbose = FALSE)
    diag_obs_mat <- unions_to_matrix(diag_obs$unions, K)
    diag_p2 <- fit_loglinear(diag_obs_mat, null_result$null_mats, K)

    diag_cs <- construct_choice_sets_m5(diag_obs$unions, men, women, J = J_CHOICE)
    diag_cl <- estimate_clogit_m5(diag_cs, K, GROUP_NAMES)

    vcat("\n  --- True alpha matrix ---\n")
    vprint(round(ALPHA, 2))

    vcat("\n  --- Observed union matrix (bilateral 50%) ---\n")
    rownames(diag_obs_mat) <- GROUP_NAMES; colnames(diag_obs_mat) <- GROUP_NAMES
    vprint(diag_obs_mat)
    vcat(sprintf("    Total unions: %d (of %d possible)\n", sum(diag_obs_mat), sum(GROUP_SIZES)))

    vcat("\n  --- Average null matrix ---\n")
    tmp_null <- round(diag_p2$avg_null, 1)
    rownames(tmp_null) <- GROUP_NAMES; colnames(tmp_null) <- GROUP_NAMES
    vprint(tmp_null)

    vcat("\n  --- Crude log(obs/null) ---\n")
    tmp_crude <- round(diag_p2$crude_ratio, 2)
    rownames(tmp_crude) <- GROUP_NAMES; colnames(tmp_crude) <- GROUP_NAMES
    vprint(tmp_crude)

    vcat("\n  --- Log-linear three-way interaction (k=1, observed) ---\n")
    tmp_3way <- round(diag_p2$three_way, 3)
    rownames(tmp_3way) <- GROUP_NAMES; colnames(tmp_3way) <- GROUP_NAMES
    vprint(tmp_3way)

    vcat("\n  --- Clogit K×K coefficient matrix ---\n")
    tmp_cl <- round(diag_cl$coef_matrix, 2)
    rownames(tmp_cl) <- GROUP_NAMES; colnames(tmp_cl) <- GROUP_NAMES
    vprint(tmp_cl)
    vcat(sprintf("    Clogit converged: %d/%d groups\n",
                 diag_cl$n_converged, K))

    vcat("\n  --- Symmetric alpha: (alpha + alpha')/2 ---\n")
    vprint(round((ALPHA + t(ALPHA))/2, 2))

    # Recovery comparison against BOTH targets
    vcat("\n  --- Recovery: all methods vs alpha_gh AND vs symmetric alpha ---\n")
    ALPHA_SYM <- (ALPHA + t(ALPHA)) / 2
    vcat(sprintf("  %-20s %10s %10s %10s %10s %10s %10s\n",
                "Method", "r(α_gh)", "ρ(α_gh)", "r(sym_α)", "ρ(sym_α)", "r_diag", "r_off"))
    vcat(paste(rep("-", 85), collapse = ""), "\n")

    for (nm_est in list(
      list("crude_ratio", diag_p2$crude_ratio),
      list("loglinear_3way", diag_p2$three_way),
      list("clogit_K2", diag_cl$coef_matrix)
    )) {
      nm <- nm_est[[1]]; est <- nm_est[[2]]
      if (any(is.na(est))) {
        vcat(sprintf("  %-20s %10s %10s %10s %10s %10s %10s\n",
                     nm, "NA", "NA", "NA", "NA", "NA", "NA"))
      } else {
        r_a <- cor(as.vector(est), as.vector(ALPHA), use = "complete.obs")
        rho_a <- cor(as.vector(est), as.vector(ALPHA), method = "spearman", use = "complete.obs")
        r_s <- cor(as.vector(est), as.vector(ALPHA_SYM), use = "complete.obs")
        rho_s <- cor(as.vector(est), as.vector(ALPHA_SYM), method = "spearman", use = "complete.obs")
        # Also diagonal and off-diagonal
        r_d <- cor(diag(est), diag(ALPHA), use = "complete.obs")
        mask <- !diag(rep(TRUE, K))
        r_o <- cor(est[mask], ALPHA[mask], use = "complete.obs")
        vcat(sprintf("  %-20s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
                     nm, r_a, rho_a, r_s, rho_s, r_d, r_o))
      }
    }

    vcat("\n  --- KEY CHECK 1: Does 3-way differ from crude ratio? ---\n")
    diff <- diag_p2$three_way - diag_p2$crude_ratio
    vcat(sprintf("    Max |diff|:  %.4f\n", max(abs(diff), na.rm = TRUE)))
    vcat(sprintf("    Mean |diff|: %.4f\n", mean(abs(diff), na.rm = TRUE)))
    vcat(sprintf("    Crude ratio range:   [%.2f, %.2f]\n",
                 min(diag_p2$crude_ratio, na.rm=TRUE), max(diag_p2$crude_ratio, na.rm=TRUE)))
    vcat(sprintf("    Three-way range:     [%.3f, %.3f]\n",
                 min(diag_p2$three_way, na.rm=TRUE), max(diag_p2$three_way, na.rm=TRUE)))
    if (max(abs(diff), na.rm = TRUE) < 0.01) {
      vcat("    *** WARNING: Three-way ≈ crude ratio — margin conditioning NOT working! ***\n")
    } else {
      vcat("    OK: Three-way differs from crude ratio.\n")
    }

    vcat("\n  --- KEY CHECK 2: loglinear() param extraction ---\n")
    vcat(sprintf("    Converged: %s\n", diag_p2$loglinear_converged))
    vcat(sprintf("    N scenarios in table: %d (expected: %d)\n",
                 diag_p2$n_scenarios, N_NULL_SIMS + 1))
    # Check if param extraction worked or fell through to manual
    # Re-run loglinear to inspect param structure
    test_tab <- array(0, dim = c(K, K, N_NULL_SIMS + 1))
    test_tab[, , 1] <- diag_obs_mat
    for (s in 1:N_NULL_SIMS) test_tab[, , s+1] <- null_result$null_mats[[s]]
    test_tab[test_tab == 0] <- 0.5
    test_fit <- tryCatch(
      loglin(test_tab, margin = list(c(1,2,3)), fit = TRUE, param = TRUE, print = FALSE),
      error = function(e) NULL
    )
    if (!is.null(test_fit) && !is.null(test_fit$param)) {
      vcat(sprintf("    loglin param names: %s\n", paste(names(test_fit$param), collapse = ", ")))
      tw <- test_fit$param[["1.2.3"]]
      if (!is.null(tw)) {
        vcat(sprintf("    Three-way param dim: %s\n", paste(dim(tw), collapse = " x ")))
        vcat("    Three-way param[,,1] (first 3 rows, 3 cols):\n")
        vprint(round(tw[1:min(3,K), 1:min(3,K), 1], 4))
      } else {
        vcat("    WARNING: param '1.2.3' not found. Names available:\n")
        vcat(sprintf("    %s\n", paste(names(test_fit$param), collapse = ", ")))
        # Show last param which should be the highest-order
        last_param <- test_fit$param[[length(test_fit$param)]]
        vcat(sprintf("    Last param dim: %s\n",
                     if(is.null(dim(last_param))) "scalar" else paste(dim(last_param), collapse=" x ")))
      }
    } else {
      vcat("    loglin() failed or param=NULL — using manual fallback extraction.\n")
    }

    vcat("\n  --- KEY CHECK 3: Clogit diagonal vs alpha diagonal ---\n")
    if (!any(is.na(diag_cl$coef_matrix))) {
      for (g in 1:K) {
        vcat(sprintf("    %s: alpha=%.2f  clogit=%.2f  crude=%.2f  3way=%.3f\n",
                     GROUP_NAMES[g],
                     ALPHA[g,g],
                     diag_cl$coef_matrix[g,g],
                     diag_p2$crude_ratio[g,g],
                     diag_p2$three_way[g,g]))
      }
    }

    vcat("\n  --- KEY CHECK 4: A specific asymmetric pair ---\n")
    vcat("    Nordic→MENA: true α = -0.6, MENA→Nordic: true α = -0.8\n")
    g_n <- 1; g_m <- 4  # Nordic=1, MENA=4 (in K=5 config)
    vcat(sprintf("    Crude:  Nord→MENA=%.2f  MENA→Nord=%.2f\n",
                 diag_p2$crude_ratio[g_n,g_m], diag_p2$crude_ratio[g_m,g_n]))
    vcat(sprintf("    3-way:  Nord→MENA=%.3f  MENA→Nord=%.3f\n",
                 diag_p2$three_way[g_n,g_m], diag_p2$three_way[g_m,g_n]))
    if (!any(is.na(diag_cl$coef_matrix))) {
      vcat(sprintf("    Clogit: Nord→MENA=%.2f  MENA→Nord=%.2f\n",
                   diag_cl$coef_matrix[g_n,g_m], diag_cl$coef_matrix[g_m,g_n]))
    }

    vcat("\n  ======== END VALIDATION DIAGNOSTICS ========\n")
    vcat(sprintf("  Validation log saved to: %s\n\n", val_log))
    close(val_con)
  }

  # =========================================================================
  # ONE-SIDED SANITY CHECK (first rep only)
  # =========================================================================
  if (rep == 1) {
    cat("  One-sided sanity check...\n")
    t0 <- Sys.time()

    os_result <- onesided_gravity_match(men, women, ALPHA, DELTA)
    os_mat <- unions_to_matrix(os_result$unions, K)
    p2_os <- fit_loglinear(os_mat, null_result$null_mats, K)

    rec_os_crude <- compute_recovery(p2_os$crude_ratio, ALPHA, K, GROUP_NAMES)
    rec_os_3way  <- compute_recovery(p2_os$three_way, ALPHA, K, GROUP_NAMES)

    cs_os <- construct_choice_sets_m5(os_result$unions, men, women, J = J_CHOICE)
    clogit_os <- estimate_clogit_m5(cs_os, K, GROUP_NAMES)
    if (clogit_os$n_converged >= K - 1) {
      rec_os_clogit <- compute_recovery(clogit_os$coef_matrix, ALPHA, K, GROUP_NAMES)
    } else {
      rec_os_clogit <- list(r_full = NA, rho_full = NA, r_diag = NA, rho_diag = NA,
                             r_off = NA, rho_off = NA, r_sym = NA, r_min = NA,
                             rmse_full = NA, rmse_diag = NA)
    }

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    for (method_info in list(
      list(name = "crude_ratio",    rec = rec_os_crude),
      list(name = "loglinear_3way", rec = rec_os_3way),
      list(name = "clogit_K2",     rec = rec_os_clogit)
    )) {
      row <- data.frame(
        rep = 1, condition = "onesided", acceptance_rate = NA,
        method = method_info$name,
        r_full = method_info$rec$r_full, rho_full = method_info$rec$rho_full,
        r_diag = method_info$rec$r_diag, rho_diag = method_info$rec$rho_diag,
        r_off = method_info$rec$r_off, rho_off = method_info$rec$rho_off,
        r_sym = method_info$rec$r_sym, r_min = method_info$rec$r_min,
        rmse_full = method_info$rec$rmse_full, rmse_diag = method_info$rec$rmse_diag,
        pairing_rate = os_result$pairing_rate,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
      write.table(row, results_file, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    }

    cat(sprintf("    One-sided sanity check done (%.0f sec).\n", elapsed))
    cat("\n    SANITY CHECK RESULTS:\n")
    print_recovery_comparison(
      list("Crude ratio" = rec_os_crude,
           "Log-linear 3-way" = rec_os_3way,
           "Clogit K²" = rec_os_clogit),
      ALPHA, GROUP_NAMES
    )
  }

  # Memory cleanup
  gc()

  cat(sprintf("  Rep %d total: %.0f sec\n", rep,
              as.numeric(difftime(Sys.time(), t0_rep, units = "secs"))))
}


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("FINAL SUMMARY\n")
cat(rep("=", 70), "\n\n", sep = "")

# Read full results
all_results <- read.csv(results_file)

# Aggregate by condition × method
cat("Mean recovery metrics across MC reps:\n\n")

for (cond in unique(all_results$condition)) {
  for (ar in unique(all_results$acceptance_rate[all_results$condition == cond])) {
    subset <- all_results[all_results$condition == cond &
                          (is.na(ar) | all_results$acceptance_rate == ar), ]
    if (nrow(subset) == 0) next

    label <- if (cond == "bilateral") {
      sprintf("Bilateral (%.0f%% acceptance)", ar * 100)
    } else {
      "One-sided (sanity check)"
    }
    cat(sprintf("\n--- %s ---\n", label))
    cat(sprintf("%-20s %8s %8s %8s %8s %8s %8s %8s\n",
                "Method", "r_full", "ρ_full", "r_diag", "r_off", "ρ_off", "r_sym", "RMSE"))
    cat(paste(rep("-", 80), collapse = ""), "\n")

    for (m in unique(subset$method)) {
      ms <- subset[subset$method == m, ]
      cat(sprintf("%-20s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                  m,
                  mean(ms$r_full, na.rm = TRUE),
                  mean(ms$rho_full, na.rm = TRUE),
                  mean(ms$r_diag, na.rm = TRUE),
                  mean(ms$r_off, na.rm = TRUE),
                  mean(ms$rho_off, na.rm = TRUE),
                  mean(ms$r_sym, na.rm = TRUE),
                  mean(ms$rmse_full, na.rm = TRUE)))
    }
  }
}

cat(sprintf("\n\nResults saved to: %s\n", results_file))
cat("Done.\n")
