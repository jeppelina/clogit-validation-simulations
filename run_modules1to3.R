# =============================================================================
# Production run: Modules 1, 1b, 2, 3 — Crash-safe sequential execution
# =============================================================================
#
# INCREMENTAL SAVE DESIGN: Results are appended to CSV after every single
# condition. If R crashes, you keep everything that completed.
#
# Uses the same pattern as run_module4_production.R:
#   - Sequential execution only (N_CORES = 1)
#   - Incremental CSV saves after each condition
#   - Aggressive gc() to prevent memory accumulation
#
# Usage:
#   Rscript run_modules1to3.R                         # all modules
#   Rscript run_modules1to3.R module1                 # Module 1 only
#   Rscript run_modules1to3.R module1b                # Module 1b only
#   Rscript run_modules1to3.R all quick               # all, reduced reps
#
# Output:
#   results/module1_results.csv
#   results/module1b_results.csv
#
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- ("quick" %in% args)
run_module <- if (length(setdiff(args, "quick")) > 0) setdiff(args, "quick")[1] else "all"

source("utils.R")
library(nnet)  # for module 1b mlogit

# =============================================================================
# PARAMETERS
# =============================================================================

N_CORES <- 1  # sequential only — mclapply corrupts memory on macOS with clogit

if (quick_mode) {
  R_SIMS   <- 20
  N_AGENTS <- 500
  cat("=== QUICK MODE: R=20, N=500 ===\n\n")
} else {
  R_SIMS   <- 200
  N_AGENTS <- 1000
  cat(sprintf("=== PRODUCTION: R=%d, N=%d ===\n\n", R_SIMS, N_AGENTS))
}

J_ALTS <- 30

# Shared DGP parameters
BETA_EDU        <- -1.5   # homophily coefficient (Module 1, 1b)
BETA_SAME_GROUP <- 0.8    # endogamy coefficient
BETA_STATUS     <- 0.5    # status max coefficient (Module 1b)

dir.create("results", showWarnings = FALSE)


# =============================================================================
# CRASH-SAFE INFRASTRUCTURE
# =============================================================================

append_results <- function(df, file) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  write_header <- !file.exists(file)
  write.table(df, file = file, sep = ",", row.names = FALSE,
              col.names = write_header, append = !write_header)
  cat(sprintf("  >> Saved %d rows to %s\n", nrow(df), file))
}

run_mc_safe <- function(sim_fn, R, seed = 42, ...) {
  results <- vector("list", R)
  for (r in 1:R) {
    set.seed(seed + r)
    results[[r]] <- tryCatch(sim_fn(...), error = function(e) NULL)
    if (r %% 50 == 0) {
      cat(sprintf("    [%d/%d]\n", r, R))
      gc(verbose = FALSE)
    }
  }
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) return(NULL)
  cat(sprintf("  %d / %d valid\n", length(results), R))
  results
}

# Generic: run MC, assess recovery, save, free memory
run_condition <- function(module, condition_label, sim_fn, true_params,
                          model_name, extra_cols = list(),
                          results_file, ...) {
  cat(sprintf("\n--- [%s] %s ---\n", module, condition_label))
  t0 <- proc.time()["elapsed"]

  mc <- tryCatch(
    run_mc_safe(sim_fn, R = R_SIMS, seed = 42, ...),
    error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
  )

  elapsed <- proc.time()["elapsed"] - t0

  if (is.null(mc) || length(mc) < 5) {
    cat(sprintf("  SKIPPED (%s valid, %.1fs)\n",
                if (is.null(mc)) "0" else as.character(length(mc)), elapsed))
    gc(verbose = FALSE)
    return(invisible(NULL))
  }

  # Filter to runs that have the requested model
  mc_filt <- mc[sapply(mc, function(x) is.list(x) && !is.null(x[[model_name]]))]
  if (length(mc_filt) < 5) {
    cat(sprintf("  WARNING: only %d runs have '%s', skipping\n",
                length(mc_filt), model_name))
    gc(verbose = FALSE)
    return(invisible(NULL))
  }

  rec <- assess_recovery(mc_filt, true_params, model_name)

  rows <- data.frame(
    module     = module,
    condition  = condition_label,
    model      = model_name,
    parameter  = rec$parameter,
    true_value = rec$true_value,
    mean_est   = rec$mean_estimate,
    median_est = rec$median_estimate,
    sd_est     = rec$sd_estimate,
    bias       = rec$bias,
    rel_bias_pct = rec$rel_bias_pct,
    rmse       = rec$rmse,
    coverage_95 = rec$coverage_95,
    n_valid    = rec$n_valid,
    n_runs     = R_SIMS,
    elapsed_min = round(elapsed / 60, 2),
    stringsAsFactors = FALSE
  )

  for (nm in names(extra_cols)) {
    rows[[nm]] <- extra_cols[[nm]]
  }

  for (j in 1:nrow(rows)) {
    cat(sprintf("  %s: bias=%+.4f (%+.1f%%)  rmse=%.4f  cov95=%.3f\n",
                rows$parameter[j], rows$bias[j], rows$rel_bias_pct[j],
                rows$rmse[j], rows$coverage_95[j]))
  }
  cat(sprintf("  (%.1f min)\n", elapsed / 60))

  append_results(rows, results_file)

  rm(mc, mc_filt, rec)
  gc(verbose = FALSE)

  invisible(rows)
}


# =============================================================================
# MODULE 1: Group Size and Uneven Trait Distributions
# =============================================================================

run_module1_production <- function() {
  cat("\n\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("MODULE 1: Group Size and Uneven Trait Distributions\n")
  cat("=", rep("=", 70), "\n", sep = "")

  RESULTS_FILE <- "results/module1_results.csv"
  if (file.exists(RESULTS_FILE)) file.remove(RESULTS_FILE)

  experiments <- list(
    balanced = list(
      name = "Balanced (33/33/33)",
      group_proportions = c(1/3, 1/3, 1/3),
      edu_means = c(12, 12, 12), edu_sds = c(2, 2, 2)
    ),
    moderate = list(
      name = "Moderate (50/30/20)",
      group_proportions = c(0.50, 0.30, 0.20),
      edu_means = c(12, 10, 14), edu_sds = c(2, 3, 1.5)
    ),
    unbalanced = list(
      name = "Unbalanced (70/20/10)",
      group_proportions = c(0.70, 0.20, 0.10),
      edu_means = c(12, 10, 14), edu_sds = c(2, 3, 1.5)
    ),
    extreme = list(
      name = "Extreme (85/10/5)",
      group_proportions = c(0.85, 0.10, 0.05),
      edu_means = c(12, 9, 15), edu_sds = c(2, 3, 1)
    )
  )

  # Sim function for Module 1 — oracle (one-sided)
  m1_sim <- function(group_proportions, edu_means, edu_sds) {
    agents <- generate_agents(
      N = N_AGENTS,
      n_groups = length(group_proportions),
      group_proportions = group_proportions,
      edu_means = edu_means,
      edu_sds = edu_sds
    )
    gen_beta <- c(BETA_EDU, BETA_SAME_GROUP)
    gen_vars <- c("edu_diff", "same_group")
    choice_data <- construct_choice_sets(agents, J = J_ALTS)
    choice_data <- compute_match_variables(choice_data)
    choice_data <- simulate_choices(choice_data, beta = gen_beta, var_names = gen_vars)

    res_A <- estimate_clogit(choice_data, "edu_diff + same_group")
    choice_data$edu_diff_x_same <- choice_data$edu_diff * choice_data$same_group
    res_B <- estimate_clogit(choice_data, "edu_diff + same_group + edu_diff_x_same")

    list(model_A = res_A, model_B = res_B)
  }

  # Sim function for Module 1 — bilateral matching
  m1_sim_bilateral <- function(group_proportions, edu_means, edu_sds,
                                bilateral_threshold = -0.5,
                                bilateral_rounds = 200) {
    agents <- generate_agents(
      N = N_AGENTS,
      n_groups = length(group_proportions),
      group_proportions = group_proportions,
      edu_means = edu_means,
      edu_sds = edu_sds
    )
    gen_beta <- c(BETA_EDU, BETA_SAME_GROUP)
    gen_vars <- c("edu_diff", "same_group")

    bm <- bilateral_matching(agents, gen_beta, gen_vars,
                             n_rounds = bilateral_rounds,
                             threshold = bilateral_threshold)
    if (nrow(bm$matches) < 20) return(NULL)

    choice_data <- construct_choice_sets_bilateral(agents, bm$matches, J = J_ALTS)
    choice_data <- compute_match_variables(choice_data)

    res_A <- estimate_clogit(choice_data, "edu_diff + same_group")

    list(model_A = res_A)
  }

  true_A <- c(edu_diff = BETA_EDU, same_group = BETA_SAME_GROUP)
  true_B <- c(edu_diff = BETA_EDU, same_group = BETA_SAME_GROUP, edu_diff_x_same = 0)

  for (exp_name in names(experiments)) {
    exp <- experiments[[exp_name]]

    # --- Oracle ---
    run_condition(
      module = "M1", condition_label = exp$name,
      sim_fn = m1_sim, true_params = true_A, model_name = "model_A",
      extra_cols = list(groups = paste(round(exp$group_proportions * 100), collapse = "/"),
                        mechanism = "oracle"),
      results_file = RESULTS_FILE,
      group_proportions = exp$group_proportions,
      edu_means = exp$edu_means,
      edu_sds = exp$edu_sds
    )

    run_condition(
      module = "M1", condition_label = paste(exp$name, "(interaction model)"),
      sim_fn = m1_sim, true_params = true_B, model_name = "model_B",
      extra_cols = list(groups = paste(round(exp$group_proportions * 100), collapse = "/"),
                        mechanism = "oracle"),
      results_file = RESULTS_FILE,
      group_proportions = exp$group_proportions,
      edu_means = exp$edu_means,
      edu_sds = exp$edu_sds
    )

    # --- Bilateral (calibrated at 20% acceptance rate) ---
    # Calibrate threshold for this specific group composition
    cal_agents <- generate_agents(
      N = N_AGENTS, n_groups = length(exp$group_proportions),
      group_proportions = exp$group_proportions,
      edu_means = exp$edu_means, edu_sds = exp$edu_sds
    )
    cal <- calibrate_threshold(cal_agents,
                               c(BETA_EDU, BETA_SAME_GROUP),
                               c("edu_diff", "same_group"),
                               target_acceptance_rate = 0.20)
    bi_rounds <- min(ceiling(20 / 0.20^2), 500)
    rm(cal_agents); gc(verbose = FALSE)

    cat(sprintf("  Bilateral threshold for %s: %.3f (target 20%% accept)\n",
                exp$name, cal$threshold))

    run_condition(
      module = "M1", condition_label = paste(exp$name, "(bilateral 20%)"),
      sim_fn = m1_sim_bilateral, true_params = true_A, model_name = "model_A",
      extra_cols = list(groups = paste(round(exp$group_proportions * 100), collapse = "/"),
                        mechanism = "bilateral_20pct",
                        threshold = cal$threshold),
      results_file = RESULTS_FILE,
      group_proportions = exp$group_proportions,
      edu_means = exp$edu_means,
      edu_sds = exp$edu_sds,
      bilateral_threshold = cal$threshold,
      bilateral_rounds = bi_rounds
    )
  }

  cat("\n  Module 1 complete.\n")
}


# =============================================================================
# MODULE 1b: Identification — Homophily vs Status Maximisation (Part A)
# =============================================================================

run_module1b_production <- function() {
  cat("\n\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("MODULE 1b: Identification — Homophily vs Status Maximisation\n")
  cat("=", rep("=", 70), "\n", sep = "")

  RESULTS_FILE <- "results/module1b_results.csv"
  if (file.exists(RESULTS_FILE)) file.remove(RESULTS_FILE)

  GROUP_PROPS <- c(0.70, 0.20, 0.10)
  EDU_MEANS   <- c(12, 10, 14)
  EDU_SDS     <- c(2, 3, 1.5)

  dgp_types <- c("homophily", "status_max", "mixed")
  spec_info <- list(
    spec_H = list(formula = "edu_diff + same_group",
                  true_homophily = c(edu_diff = BETA_EDU, same_group = BETA_SAME_GROUP),
                  true_status    = c(edu_diff = 0, same_group = BETA_SAME_GROUP),
                  true_mixed     = c(edu_diff = BETA_EDU, same_group = BETA_SAME_GROUP)),
    spec_S = list(formula = "alt_edu + same_group",
                  true_homophily = c(alt_edu = 0, same_group = BETA_SAME_GROUP),
                  true_status    = c(alt_edu = BETA_STATUS, same_group = BETA_SAME_GROUP),
                  true_mixed     = c(alt_edu = BETA_STATUS, same_group = BETA_SAME_GROUP)),
    spec_M = list(formula = "edu_diff + alt_edu + same_group",
                  true_homophily = c(edu_diff = BETA_EDU, alt_edu = 0, same_group = BETA_SAME_GROUP),
                  true_status    = c(edu_diff = 0, alt_edu = BETA_STATUS, same_group = BETA_SAME_GROUP),
                  true_mixed     = c(edu_diff = BETA_EDU, alt_edu = BETA_STATUS, same_group = BETA_SAME_GROUP))
  )

  # Sim function using source file's logic (oracle, one-sided)
  m1b_sim <- function(dgp_type) {
    agents <- generate_agents(
      N = N_AGENTS,
      n_groups = length(GROUP_PROPS),
      group_proportions = GROUP_PROPS,
      edu_means = EDU_MEANS,
      edu_sds = EDU_SDS
    )

    if (dgp_type == "homophily") {
      gen_beta <- c(BETA_EDU, BETA_SAME_GROUP)
      gen_vars <- c("edu_diff", "same_group")
    } else if (dgp_type == "status_max") {
      gen_beta <- c(BETA_STATUS, BETA_SAME_GROUP)
      gen_vars <- c("alt_edu", "same_group")
    } else {
      gen_beta <- c(BETA_EDU, BETA_STATUS, BETA_SAME_GROUP)
      gen_vars <- c("edu_diff", "alt_edu", "same_group")
    }

    choice_data <- construct_choice_sets(agents, J = J_ALTS)
    choice_data <- compute_match_variables(choice_data)
    choice_data$alt_edu <- choice_data$alt_education
    choice_data$ego_x_alt_edu <- choice_data$ego_education * choice_data$alt_education
    choice_data <- simulate_choices(choice_data, beta = gen_beta, var_names = gen_vars)

    results <- list()
    for (spec_name in names(spec_info)) {
      results[[spec_name]] <- tryCatch(
        estimate_clogit(choice_data, spec_info[[spec_name]]$formula),
        error = function(e) NULL
      )
    }
    results
  }

  for (dgp in dgp_types) {
    for (spec_name in names(spec_info)) {
      # Map dgp name to spec_info key (status_max → true_status, not true_status_max)
      dgp_key <- ifelse(dgp == "status_max", "status", dgp)
      true_key <- paste0("true_", dgp_key)
      true_params <- spec_info[[spec_name]][[true_key]]
      if (is.null(true_params)) next

      label <- sprintf("DGP=%s, Spec=%s", dgp, spec_name)

      run_condition(
        module = "M1b", condition_label = label,
        sim_fn = m1b_sim, true_params = true_params,
        model_name = spec_name,
        extra_cols = list(dgp = dgp, spec = spec_name),
        results_file = RESULTS_FILE,
        dgp_type = dgp
      )
    }
  }

  cat("\n  Module 1b complete.\n")
}


# =============================================================================
# MODULE 2: Spatial Opportunity vs Homophily
# =============================================================================

run_module2_production <- function() {
  cat("\n\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("MODULE 2: Spatial Opportunity vs Homophily\n")
  cat("=", rep("=", 70), "\n", sep = "")

  RESULTS_FILE <- "results/module2_results.csv"
  if (file.exists(RESULTS_FILE)) file.remove(RESULTS_FILE)

  # We run Module 2 by sourcing its sim function.
  # But to keep it crash-safe, we inline a simplified version here.

  GRID_SIZE <- 100
  BETA_HOMOPHILY_M2  <- -1.5
  BETA_PROPINQUITY   <- -2.0
  BETA_SAME_GROUP_M2 <- 0.8

  # Source the module to get its sim function
  # (We can't directly source because it would execute the whole script)
  # Instead, we define the experiments and call the infrastructure

  cat("  NOTE: Module 2 uses the existing module2_spatial.R script directly.\n")
  cat("  It already has R_SIMS=200 and sequential execution.\n")
  cat("  Run with: Rscript module2_spatial.R\n")
  cat("  Module 2 skipped in this runner — use the standalone script.\n")
}


# =============================================================================
# MODULE 3: Mediation via Sequential Opportunity Variables
# =============================================================================

run_module3_production <- function() {
  cat("\n\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("MODULE 3: Mediation via Sequential Opportunity Variables\n")
  cat("=", rep("=", 70), "\n", sep = "")

  RESULTS_FILE <- "results/module3_results.csv"
  if (file.exists(RESULTS_FILE)) file.remove(RESULTS_FILE)

  # Module 3 parameters
  BETA_EDU_M3   <- -1.0
  BETA_WP       <- 1.5
  BETA_OPENNESS <- -0.8
  EDU_TO_WP     <- 2.0
  OPENNESS_TO_WP <- 1.0
  N_WP          <- 20

  true_M1 <- c(edu_diff = BETA_EDU_M3, same_wp = BETA_WP)

  # We need Module 3's sim function. Since it's complex and already in
  # module3_mediation.R, the safest approach is to run that script directly
  # with the updated R_SIMS=200 and oracle-only settings.

  cat("  NOTE: Module 3 is complex (Parts I-III, 8 assumption sweeps).\n")
  cat("  It already has R_SIMS=200, N_CORES=1, and oracle-only matching.\n")
  cat("  Run with: Rscript module3_mediation.R\n")
  cat("  \n")
  cat("  IMPORTANT: Module 3 saves a large .rds file. To prevent OOM,\n")
  cat("  consider running parts separately or adding incremental saves\n")
  cat("  to the script. Monitor memory usage during execution.\n")
  cat("  Module 3 skipped in this runner — use the standalone script.\n")
}


# =============================================================================
# MAIN
# =============================================================================

cat(sprintf("Run started at %s\n", Sys.time()))
cat(sprintf("Module selection: %s\n\n", run_module))

if (run_module %in% c("all", "module1")) run_module1_production()
if (run_module %in% c("all", "module1b")) run_module1b_production()
if (run_module %in% c("all", "module2")) run_module2_production()
if (run_module %in% c("all", "module3")) run_module3_production()

cat(sprintf("\n\nAll done at %s\n", Sys.time()))
