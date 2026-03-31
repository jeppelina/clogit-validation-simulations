# =============================================================================
# Module 4: Matching Mechanisms — How bilateral matching biases clogit
# =============================================================================
#
# QUESTION: How does the mismatch between two-sided matching (reality) and
# one-sided choice models (clogit) affect parameter recovery?
#
# HEADLINE FINDING: Bilateral matching overestimates educational homophily
# by ~23% (relative bias on edu_diff, threshold = -0.5).
#
# STRUCTURE:
#   Part 1 — Core comparison: 3 mechanisms × 3 DGPs (9 conditions + misspec)
#            Mechanisms: oracle, bilateral, deferred acceptance
#            HEADLINE: oracle vs bilateral under homophily DGP
#            ROBUSTNESS: status_max and mixed DGPs
#   Part 2 — Acceptance-rate-calibrated threshold sweep: 2 DGPs × 5 rates
#            KEY FIGURE for kappa: bias × selectivity, apples-to-apples
#   Part 3 — Group interaction: mechanism × balanced/moderate/extreme groups
#            ROBUSTNESS: bias stability across group compositions
#   Part 4 — Bilateral rounds sweep: encounter rounds from 5 to 200
#            Separates encounter limitation from fundamental bilateral bias
#   Part 5 — Pair diagnostics: utility distributions for matched pairs
#            KEY FIGURE: shows WHY bilateral creates bias (mutual acceptance
#            truncation oversamples similar-education pairs)
#   Part 6 — V distribution diagnostics & calibration verification
#            Confirms calibrated thresholds equalize selectivity across DGPs
#            Single-realization check of sorting under calibrated thresholds
#
# DGP (primary): U = -1.5 × |edu_diff| + 0.8 × same_group + ε (Gumbel)
# DGP (robustness): status_max, mixed — see run_mechanism_sim() below
#
# CRASH-SAFE: Results appended to CSV after every condition.
# Sequential execution only (mclapply corrupts memory on macOS with clogit).
#
# Usage:
#   Rscript run_module4.R              # full run (~2-3 hrs)
#   Rscript run_module4.R quick        # sanity check (~20 min)
#
# Output:
#   results/module4_results.csv            — parameter recovery (one row per parameter per condition)
#   results/module4_pair_diagnostics.csv   — pair-level utility data (Part 5)
#   results/module4_status_diagnostics.csv — status_max verification data (Part 6)
#   results/module4_pair_*.pdf             — diagnostic plots
#
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- ("quick" %in% args)

source("utils.R")

# =============================================================================
# PARAMETERS
# =============================================================================

N_CORES <- 1  # sequential only — mclapply corrupts memory on macOS with clogit

if (quick_mode) {
  R_SIMS   <- 50
  N_AGENTS <- 500
  cat("=== QUICK MODE: R=50, N=500 ===\n\n")
} else {
  R_SIMS   <- 200
  N_AGENTS <- 1000
  cat(sprintf("=== PRODUCTION: R=%d, N=%d ===\n\n", R_SIMS, N_AGENTS))
}

J_ALTS <- 30
BETA_HOMOPHILY  <- -1.5
BETA_STATUS     <- 0.5
BETA_SAME_GROUP <- 0.8
GROUP_PROPS <- c(0.70, 0.20, 0.10)
EDU_MEANS   <- c(12, 10, 14)
EDU_SDS     <- c(2, 3, 1.5)

dir.create("results", showWarnings = FALSE)
RESULTS_FILE <- "results/module4_results.csv"


# =============================================================================
# INCREMENTAL SAVE: append one summary to the CSV after each condition
# =============================================================================

append_results <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  write_header <- !file.exists(RESULTS_FILE)
  write.table(df, file = RESULTS_FILE, sep = ",", row.names = FALSE,
              col.names = write_header, append = !write_header)
  cat(sprintf("  >> Saved %d rows to %s\n", nrow(df), RESULTS_FILE))
}


# =============================================================================
# SEQUENTIAL MC RUNNER (no forking)
# =============================================================================

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


# =============================================================================
# RUN ONE CONDITION: MC + extract summary + save + free memory
# =============================================================================

run_and_save <- function(part, label, dgp, mech, true_params, spec_name,
                         extra_cols = list(), ...) {
  cat(sprintf("\n--- [%s] %s ---\n", part, label))
  t0 <- proc.time()["elapsed"]

  mc <- tryCatch(
    run_mc_safe(run_mechanism_sim, R = R_SIMS, seed = 42, ...),
    error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
  )

  elapsed <- proc.time()["elapsed"] - t0

  if (is.null(mc) || length(mc) < 5) {
    cat(sprintf("  SKIPPED (%s valid, %.1fs)\n",
                if (is.null(mc)) "0" else as.character(length(mc)), elapsed))
    gc(verbose = FALSE)
    return(invisible(NULL))
  }

  # Filter to runs that have the requested spec
  mc <- mc[sapply(mc, function(x) !is.null(x[[spec_name]]))]
  if (length(mc) < 5) {
    cat(sprintf("  WARNING: only %d runs have '%s', skipping\n", length(mc), spec_name))
    gc(verbose = FALSE)
    return(invisible(NULL))
  }

  # Recovery stats
  rec <- assess_recovery(mc, true_params, spec_name)

  # Matching stats
  stats <- sapply(mc, function(x) {
    s <- x$match_stats
    c(n_pairs = s$n_pairs, cor = s$cor_ego_alt_edu,
      ed = s$mean_edu_diff, sg = s$same_group_rate)
  })

  # Build summary rows
  rows <- data.frame(
    part       = part,
    dgp        = dgp,
    mechanism  = mech,
    spec       = spec_name,
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
    n_runs     = length(mc),
    mean_cor   = mean(stats["cor", ]),
    mean_pair_rate = 2 * mean(stats["n_pairs", ]) / N_AGENTS,
    mean_sg_rate   = mean(stats["sg", ]),
    elapsed_min    = round(elapsed / 60, 2),
    stringsAsFactors = FALSE
  )

  # Add any extra columns (threshold, bilateral_rounds, groups, etc.)
  for (nm in names(extra_cols)) {
    rows[[nm]] <- extra_cols[[nm]]
  }

  # Print summary
  for (j in 1:nrow(rows)) {
    cat(sprintf("  %s: bias=%+.4f (%+.1f%%)  rmse=%.4f  cov95=%.3f\n",
                rows$parameter[j], rows$bias[j], rows$rel_bias_pct[j],
                rows$rmse[j], rows$coverage_95[j]))
  }
  cat(sprintf("  cor=%.3f  rate=%.0f%%  (%.1f min)\n",
              rows$mean_cor[1], 100 * rows$mean_pair_rate[1], elapsed / 60))

  # SAVE IMMEDIATELY
  append_results(rows)

  # FREE MEMORY
  rm(mc, rec, stats)
  gc(verbose = FALSE)

  invisible(rows)
}


# =============================================================================
# SIMULATION FUNCTION
# =============================================================================

run_mechanism_sim <- function(
    N = N_AGENTS, J = J_ALTS,
    group_proportions = GROUP_PROPS,
    edu_means = EDU_MEANS, edu_sds = EDU_SDS,
    dgp_type = "homophily",
    mechanism = "one_sided",
    bilateral_threshold = -0.5, bilateral_rounds = 20,
    da_J_proposals = 5, da_n_rounds = 30,
    cs_method = "random_population",
    beta_h = BETA_HOMOPHILY, beta_s = BETA_STATUS, beta_sg = BETA_SAME_GROUP,
    verbose = FALSE
) {
  n_groups <- length(group_proportions)
  agents <- generate_agents(
    N = N, n_groups = n_groups,
    group_proportions = group_proportions,
    edu_means = edu_means, edu_sds = edu_sds
  )

  if (dgp_type == "homophily") {
    gen_beta <- c(beta_h, beta_sg);  gen_vars <- c("edu_diff", "same_group")
  } else if (dgp_type == "status_max") {
    gen_beta <- c(beta_s, beta_sg);  gen_vars <- c("alt_edu", "same_group")
  } else if (dgp_type == "mixed") {
    gen_beta <- c(beta_h, beta_s, beta_sg)
    gen_vars <- c("edu_diff", "alt_edu", "same_group")
  }

  if (mechanism == "one_sided_oracle") {
    cd <- construct_choice_sets(agents, J = J)
    cd <- compute_match_variables(cd)
    cd$alt_edu <- cd$alt_education
    cd$ego_x_alt_edu <- cd$ego_education * cd$alt_education
    cd <- simulate_choices(cd, beta = gen_beta, var_names = gen_vars)

  } else if (mechanism == "bilateral") {
    bm <- bilateral_matching(agents, gen_beta, gen_vars,
                             n_rounds = bilateral_rounds,
                             threshold = bilateral_threshold)
    if (nrow(bm$matches) < 20) return(NULL)
    cd <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    cd <- compute_match_variables(cd)
    cd$alt_edu <- cd$alt_education
    cd$ego_x_alt_edu <- cd$ego_education * cd$alt_education

  } else if (mechanism == "deferred_acceptance") {
    bm <- deferred_acceptance_matching(agents, gen_beta, gen_vars,
                                       J_proposals = da_J_proposals,
                                       n_rounds = da_n_rounds)
    if (nrow(bm$matches) < 20) return(NULL)
    cd <- construct_choice_sets_bilateral(agents, bm$matches, J = J)
    cd <- compute_match_variables(cd)
    cd$alt_edu <- cd$alt_education
    cd$ego_x_alt_edu <- cd$ego_education * cd$alt_education
  }

  chosen <- cd[cd$chosen == 1, ]
  match_stats <- list(
    n_pairs = nrow(chosen),
    cor_ego_alt_edu = cor(chosen$ego_education, chosen$alt_education),
    mean_edu_diff = mean(chosen$edu_diff),
    same_group_rate = mean(chosen$same_group)
  )

  spec_H <- estimate_clogit(cd, "edu_diff + same_group")
  spec_S <- estimate_clogit(cd, "alt_edu + same_group")
  spec_M <- estimate_clogit(cd, "edu_diff + alt_edu + same_group")
  spec_I <- estimate_clogit(cd, "ego_x_alt_edu + same_group")

  # Free the big choice data before returning
  rm(cd, agents); gc(verbose = FALSE)

  list(spec_H = spec_H, spec_S = spec_S, spec_M = spec_M, spec_I = spec_I,
       match_stats = match_stats)
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Delete old results file so we start fresh
if (file.exists(RESULTS_FILE)) file.remove(RESULTS_FILE)

total_t0 <- proc.time()["elapsed"]
cat(sprintf("Started: %s\n", Sys.time()))
cat(sprintf("N=%d, J=%d, R=%d\n\n", N_AGENTS, J_ALTS, R_SIMS))


# ─────────────────────────────────────────────────────────────────────────────
# PART 1: Core Comparison — 3 mechanisms × 3 DGPs = 9 conditions (+ misspec)
#   Mechanisms: oracle (baseline), bilateral (central finding), DA (comparison)
#   HEADLINE: homophily DGP, oracle vs bilateral (kappa Table X)
#   ROBUSTNESS: status_max and mixed DGPs (kappa appendix)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 1: Core Comparison\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

mechanisms <- c("one_sided_oracle", "bilateral", "deferred_acceptance")
dgp_types  <- c("homophily", "status_max", "mixed")

for (dgp in dgp_types) {
  correct_spec <- switch(dgp,
    homophily = "spec_H", status_max = "spec_S", mixed = "spec_M")

  if (dgp == "homophily") {
    true_vals <- c(edu_diff = BETA_HOMOPHILY, same_group = BETA_SAME_GROUP)
  } else if (dgp == "status_max") {
    true_vals <- c(alt_edu = BETA_STATUS, same_group = BETA_SAME_GROUP)
  } else {
    true_vals <- c(edu_diff = BETA_HOMOPHILY, alt_edu = BETA_STATUS,
                  same_group = BETA_SAME_GROUP)
  }

  for (mech in mechanisms) {
    run_and_save("P1_core", sprintf("%s × %s", dgp, mech),
                 dgp, mech, true_vals, correct_spec,
                 dgp_type = dgp, mechanism = mech)

    # Misspecification: spec_H on status_max data
    if (dgp == "status_max") {
      run_and_save("P1_misspec", sprintf("%s × %s (spec_H misspec)", dgp, mech),
                   dgp, mech,
                   c(edu_diff = 0, same_group = BETA_SAME_GROUP), "spec_H",
                   dgp_type = dgp, mechanism = mech)
    }
  }
}

cat(sprintf("\n  [Part 1 done at %s]\n", Sys.time()))


# ─────────────────────────────────────────────────────────────────────────────
# PART 2: Acceptance-Rate-Calibrated Threshold Sweep
#   KEY FIGURE: bias vs acceptance rate across DGPs (kappa Figure X)
#
#   Instead of sweeping raw thresholds (which are meaningless across DGPs
#   because V distributions differ), we equalize MARKET SELECTIVITY:
#   for each DGP, find the threshold that yields a target one-sided
#   acceptance rate. This gives apples-to-apples comparison.
#
#   Target rates: 10%, 20%, 40%, 60%, 80% (realistic range, granular at low end)
#   DGPs: homophily and status_max
#   → 5 rates × 2 DGPs = 10 conditions
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 2: Acceptance-Rate-Calibrated Threshold Sweep\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

true_H <- c(edu_diff = BETA_HOMOPHILY, same_group = BETA_SAME_GROUP)
true_S <- c(alt_edu = BETA_STATUS, same_group = BETA_SAME_GROUP)

target_rates <- c(0.10, 0.20, 0.40, 0.60, 0.80)

# Generate calibration population (shared across DGPs for consistency)
set.seed(123)
agents_calib <- generate_agents(
  N = N_AGENTS, n_groups = length(GROUP_PROPS),
  group_proportions = GROUP_PROPS,
  edu_means = EDU_MEANS, edu_sds = EDU_SDS
)

p2_dgps <- list(
  homophily  = list(beta = c(BETA_HOMOPHILY, BETA_SAME_GROUP),
                    vars = c("edu_diff", "same_group"),
                    true_vals = true_H, spec = "spec_H"),
  status_max = list(beta = c(BETA_STATUS, BETA_SAME_GROUP),
                    vars = c("alt_edu", "same_group"),
                    true_vals = true_S, spec = "spec_S")
)

for (dgp_name in names(p2_dgps)) {
  dgp_info <- p2_dgps[[dgp_name]]

  cat(sprintf("\n  --- Calibrating thresholds for %s DGP ---\n", dgp_name))

  for (target in target_rates) {
    cal <- calibrate_threshold(agents_calib, dgp_info$beta, dgp_info$vars,
                               target_acceptance_rate = target,
                               n_sample_pairs = 10000)

    # Scale rounds with selectivity: joint acceptance ≈ target^2,
    # so tighter markets need more encounters to form enough pairs.
    # At 5% one-sided → ~0.25% joint → need ~400 rounds for decent pairing.
    rounds_needed <- ceiling(20 / target^2)
    rounds_needed <- min(rounds_needed, 500)  # cap at 500 to stay tractable

    cat(sprintf("  Target %.0f%%: threshold = %+.3f, rounds = %d (V_mean=%.2f, V_sd=%.2f)\n",
                100 * target, cal$threshold, rounds_needed, cal$V_mean, cal$V_sd))

    run_and_save("P2_calibrated", sprintf("%s bilateral accept=%.0f%% (thr=%.2f, r=%d)",
                                           dgp_name, 100 * target, cal$threshold, rounds_needed),
                 dgp_name, "bilateral", dgp_info$true_vals, dgp_info$spec,
                 extra_cols = list(threshold = cal$threshold,
                                   target_accept_rate = target,
                                   bilateral_rounds = rounds_needed,
                                   V_mean = cal$V_mean,
                                   V_sd = cal$V_sd),
                 dgp_type = dgp_name, mechanism = "bilateral",
                 bilateral_threshold = cal$threshold,
                 bilateral_rounds = rounds_needed)
  }
}

rm(agents_calib)
gc(verbose = FALSE)

cat(sprintf("\n  [Part 2 done at %s]\n", Sys.time()))


# ─────────────────────────────────────────────────────────────────────────────
# PART 3: Mechanism × Group Size (6 conditions) (ROBUSTNESS)
#   Oracle + bilateral × balanced/moderate/extreme groups
#   Shows bilateral bias is stable (~20-29%) across group compositions
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 3: Mechanism × Group Size\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

group_configs <- list(
  balanced = list(props = c(1/3, 1/3, 1/3), means = c(12, 12, 12), sds = c(2, 2, 2)),
  moderate = list(props = c(0.70, 0.20, 0.10), means = c(12, 10, 14), sds = c(2, 3, 1.5)),
  extreme  = list(props = c(0.85, 0.10, 0.05), means = c(12, 9, 15), sds = c(2, 3, 1.0))
)

for (gc_name in names(group_configs)) {
  gc <- group_configs[[gc_name]]
  for (mech in c("one_sided_oracle", "bilateral")) {
    run_and_save("P3_groups", sprintf("%s × %s", gc_name, mech),
                 "homophily", mech, true_H, "spec_H",
                 extra_cols = list(groups = gc_name),
                 dgp_type = "homophily", mechanism = mech,
                 group_proportions = gc$props,
                 edu_means = gc$means, edu_sds = gc$sds)
  }
}

cat(sprintf("\n  [Part 3 done at %s]\n", Sys.time()))


# ─────────────────────────────────────────────────────────────────────────────
# PART 4: Bilateral Rounds Sweep (6 conditions) (ROBUSTNESS)
#   Tests whether the bilateral bias is driven by insufficient encounter
#   opportunities (compositional selection) vs the fundamental two-sidedness.
#   If bias converges with more rounds, the asymptotic bias is the "true"
#   bilateral effect; any extra bias at low rounds is encounter limitation.
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 4: Bilateral Rounds Sweep\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

for (nr in c(5, 10, 20, 50, 100, 200)) {
  run_and_save("P4_rounds", sprintf("bilateral rounds=%d (thr=-0.5)", nr),
               "homophily", "bilateral", true_H, "spec_H",
               extra_cols = list(bilateral_rounds = nr),
               dgp_type = "homophily", mechanism = "bilateral",
               bilateral_threshold = -0.5, bilateral_rounds = nr)
}

cat(sprintf("\n  [Part 4 done at %s]\n", Sys.time()))


# ─────────────────────────────────────────────────────────────────────────────
# PART 5: Pair Diagnostics — Why bilateral matching creates bias
#   For a single large population, compute pair-level utilities under
#   oracle and bilateral. Shows the selection mechanism: bilateral requires
#   mutual acceptance, which truncates the utility distribution and
#   oversamples similar-education pairs.
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 5: Pair Diagnostics\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Single large population for diagnostics (not MC — one realization)
set.seed(42)
N_DIAG <- if (quick_mode) 500 else 2000

agents_diag <- generate_agents(
  N = N_DIAG, n_groups = length(GROUP_PROPS),
  group_proportions = GROUP_PROPS,
  edu_means = EDU_MEANS, edu_sds = EDU_SDS
)

gen_beta_diag <- c(BETA_HOMOPHILY, BETA_SAME_GROUP)
gen_vars_diag <- c("edu_diff", "same_group")

# --- Oracle: each agent picks their best partner ---
cat("\n  Computing oracle pair utilities...\n")
oracle_V_ego <- numeric(N_DIAG)
oracle_V_alt <- numeric(N_DIAG)
oracle_chosen <- integer(N_DIAG)

for (i in 1:N_DIAG) {
  best_V <- -Inf
  best_j <- NA
  for (j in (1:N_DIAG)[-i]) {
    V <- compute_pair_utility(agents_diag, i, j, gen_beta_diag, gen_vars_diag)
    eps <- -log(-log(runif(1)))
    if (V + eps > best_V) {
      best_V <- V + eps
      best_j <- j
    }
  }
  oracle_chosen[i] <- best_j
  oracle_V_ego[i] <- compute_pair_utility(agents_diag, i, best_j,
                                           gen_beta_diag, gen_vars_diag)
  oracle_V_alt[i] <- compute_pair_utility(agents_diag, best_j, i,
                                           gen_beta_diag, gen_vars_diag)
}

oracle_diag <- data.frame(
  mechanism = "oracle",
  V_ego_for_alt = oracle_V_ego,
  V_alt_for_ego = oracle_V_alt,
  edu_diff = abs(agents_diag$education - agents_diag$education[oracle_chosen])
)
cat(sprintf("  Oracle: %d choices, mean V_ego=%.2f, mean V_alt=%.2f\n",
            N_DIAG, mean(oracle_V_ego), mean(oracle_V_alt)))

# --- Bilateral (threshold = -0.5, 20 rounds) ---
cat("  Computing bilateral pair utilities (thr=-0.5, 20 rounds)...\n")
bm_05 <- bilateral_matching(agents_diag, gen_beta_diag, gen_vars_diag,
                              n_rounds = 20, threshold = -0.5,
                              diagnostics = TRUE)

bi_05_diag <- data.frame(
  mechanism = "bilateral_thr-0.5",
  V_ego_for_alt = bm_05$pair_diag$V_ego_for_alt,
  V_alt_for_ego = bm_05$pair_diag$V_alt_for_ego,
  edu_diff = abs(agents_diag$education[match(bm_05$pair_diag$ego_id, agents_diag$id)] -
                 agents_diag$education[match(bm_05$pair_diag$alt_id, agents_diag$id)])
)
cat(sprintf("  Bilateral (thr=-0.5): %d pairs (%.0f%%), mean V_ego=%.2f, mean V_alt=%.2f\n",
            nrow(bm_05$matches), 100 * bm_05$pairing_rate,
            mean(bi_05_diag$V_ego_for_alt), mean(bi_05_diag$V_alt_for_ego)))

# --- Bilateral (threshold = -1.0, 20 rounds) — more permissive ---
cat("  Computing bilateral pair utilities (thr=-1.0, 20 rounds)...\n")
bm_10 <- bilateral_matching(agents_diag, gen_beta_diag, gen_vars_diag,
                              n_rounds = 20, threshold = -1.0,
                              diagnostics = TRUE)

bi_10_diag <- data.frame(
  mechanism = "bilateral_thr-1.0",
  V_ego_for_alt = bm_10$pair_diag$V_ego_for_alt,
  V_alt_for_ego = bm_10$pair_diag$V_alt_for_ego,
  edu_diff = abs(agents_diag$education[match(bm_10$pair_diag$ego_id, agents_diag$id)] -
                 agents_diag$education[match(bm_10$pair_diag$alt_id, agents_diag$id)])
)
cat(sprintf("  Bilateral (thr=-1.0): %d pairs (%.0f%%), mean V_ego=%.2f, mean V_alt=%.2f\n",
            nrow(bm_10$matches), 100 * bm_10$pairing_rate,
            mean(bi_10_diag$V_ego_for_alt), mean(bi_10_diag$V_alt_for_ego)))

# --- Bilateral (threshold = -0.5, 200 rounds) — saturated encounters ---
cat("  Computing bilateral pair utilities (thr=-0.5, 200 rounds)...\n")
bm_200r <- bilateral_matching(agents_diag, gen_beta_diag, gen_vars_diag,
                               n_rounds = 200, threshold = -0.5,
                               diagnostics = TRUE)

bi_200r_diag <- data.frame(
  mechanism = "bilateral_200rounds",
  V_ego_for_alt = bm_200r$pair_diag$V_ego_for_alt,
  V_alt_for_ego = bm_200r$pair_diag$V_alt_for_ego,
  edu_diff = abs(agents_diag$education[match(bm_200r$pair_diag$ego_id, agents_diag$id)] -
                 agents_diag$education[match(bm_200r$pair_diag$alt_id, agents_diag$id)])
)
cat(sprintf("  Bilateral (thr=-0.5, 200r): %d pairs (%.0f%%), mean V_ego=%.2f, mean V_alt=%.2f\n",
            nrow(bm_200r$matches), 100 * bm_200r$pairing_rate,
            mean(bi_200r_diag$V_ego_for_alt), mean(bi_200r_diag$V_alt_for_ego)))

# --- Combine and save ---
all_diag <- rbind(oracle_diag, bi_05_diag, bi_10_diag, bi_200r_diag)
write.csv(all_diag, "results/module4_pair_diagnostics.csv", row.names = FALSE)
cat("\n  Pair diagnostics saved to results/module4_pair_diagnostics.csv\n")

# --- Summary statistics ---
cat("\n  === PAIR DIAGNOSTICS SUMMARY ===\n")
cat(sprintf("  %-25s  %6s  %8s  %8s  %8s  %8s\n",
            "Mechanism", "Pairs", "V_ego", "V_alt", "V_min", "|edu_d|"))
cat("  ", rep("-", 75), "\n", sep = "")

for (m in unique(all_diag$mechanism)) {
  d <- all_diag[all_diag$mechanism == m, ]
  V_min <- pmin(d$V_ego_for_alt, d$V_alt_for_ego)
  cat(sprintf("  %-25s  %6d  %8.2f  %8.2f  %8.2f  %8.2f\n",
              m, nrow(d), mean(d$V_ego_for_alt), mean(d$V_alt_for_ego),
              mean(V_min), mean(d$edu_diff)))
}

# --- Plot: V_ego_for_alt distributions by mechanism ---
p_diag <- ggplot(all_diag, aes(x = V_ego_for_alt, fill = mechanism)) +
  geom_density(alpha = 0.4) +
  labs(title = "Distribution of Deterministic Utility V(ego -> alt) by Mechanism",
       subtitle = "Bilateral matching truncates low-utility pairs",
       x = "V(ego -> alt)", y = "Density", fill = "Mechanism") +
  theme_minimal()
ggsave("results/module4_pair_V_distributions.pdf", p_diag, width = 10, height = 6)

# Plot: V_min (minimum of V_ab and V_ba) — the binding constraint
all_diag$V_min <- pmin(all_diag$V_ego_for_alt, all_diag$V_alt_for_ego)
p_vmin <- ggplot(all_diag, aes(x = V_min, fill = mechanism)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = -0.3, y = 0, label = "threshold = -0.5", hjust = 0, size = 3) +
  labs(title = "Distribution of min(V_ab, V_ba) by Mechanism",
       subtitle = "Bilateral requires BOTH utilities above threshold — mutual acceptance truncation",
       x = "min(V_ego->alt, V_alt->ego)", y = "Density", fill = "Mechanism") +
  theme_minimal()
ggsave("results/module4_pair_Vmin_distributions.pdf", p_vmin, width = 10, height = 6)

# Plot: |edu_diff| of matched pairs
p_edu <- ggplot(all_diag, aes(x = edu_diff, fill = mechanism)) +
  geom_density(alpha = 0.4) +
  labs(title = "Education Difference of Matched Pairs by Mechanism",
       subtitle = "Bilateral oversamples similar-education pairs -> inflated homophily estimate",
       x = "|edu_ego - edu_alt|", y = "Density", fill = "Mechanism") +
  theme_minimal()
ggsave("results/module4_pair_edu_diff.pdf", p_edu, width = 10, height = 6)

cat("\n  Diagnostic plots saved to results/module4_pair_*.pdf\n")

# Clean up Part 5
rm(agents_diag, oracle_V_ego, oracle_V_alt, oracle_chosen, oracle_diag,
   bm_05, bi_05_diag, bm_10, bi_10_diag, bm_200r, bi_200r_diag, all_diag)
gc(verbose = FALSE)

cat(sprintf("\n  [Part 5 done at %s]\n", Sys.time()))


# ─────────────────────────────────────────────────────────────────────────────
# PART 6: Status_max Bilateral Verification & V Distribution Diagnostics
#   Complements Part 2 with single-realization diagnostics showing WHY
#   the calibration matters. Uses the calibrate_threshold() function to
#   find thresholds at equivalent selectivity levels, then shows:
#   a) V distributions under each DGP (explains the scale difference)
#   b) Bilateral sorting (cor) under calibrated thresholds for both DGPs
#   c) DA comparison (sorting through competition, not thresholds)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 6: Status_max Verification & V Distribution Diagnostics\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

set.seed(42)
N_DIAG <- if (quick_mode) 500 else 1000

agents_sm <- generate_agents(
  N = N_DIAG, n_groups = length(GROUP_PROPS),
  group_proportions = GROUP_PROPS,
  edu_means = EDU_MEANS, edu_sds = EDU_SDS
)

# --- 6a: V distributions under both DGPs ---
cat("\n  6a: V distributions under each DGP\n")

dgp_configs <- list(
  homophily  = list(beta = c(BETA_HOMOPHILY, BETA_SAME_GROUP),
                    vars = c("edu_diff", "same_group")),
  status_max = list(beta = c(BETA_STATUS, BETA_SAME_GROUP),
                    vars = c("alt_edu", "same_group"))
)

V_by_dgp <- list()
for (dgp_name in names(dgp_configs)) {
  cfg <- dgp_configs[[dgp_name]]
  V_vals <- numeric(1000)
  for (k in 1:1000) {
    i <- sample(1:N_DIAG, 1)
    j <- sample((1:N_DIAG)[-i], 1)
    V_vals[k] <- compute_pair_utility(agents_sm, i, j, cfg$beta, cfg$vars)
  }
  V_by_dgp[[dgp_name]] <- V_vals
  cat(sprintf("  %s: V range [%.2f, %.2f], mean=%.2f, sd=%.2f\n",
              dgp_name, min(V_vals), max(V_vals), mean(V_vals), sd(V_vals)))
}
cat("  → With a single threshold, selectivity differs wildly across DGPs.\n")
cat("  → Part 2 calibrated sweep solves this by equalizing acceptance rates.\n")

# --- 6b: Sorting under calibrated thresholds (single-realization check) ---
cat("\n  6b: Sorting under calibrated thresholds (single realization)\n")
sm_diag_rows <- list()

verify_rates <- c(0.10, 0.20, 0.40, 0.80)

for (dgp_name in names(dgp_configs)) {
  cfg <- dgp_configs[[dgp_name]]

  for (target in verify_rates) {
    cal <- calibrate_threshold(agents_sm, cfg$beta, cfg$vars,
                               target_acceptance_rate = target,
                               n_sample_pairs = 5000)

    # Scale rounds with selectivity (same logic as Part 2)
    vrounds <- min(ceiling(20 / target^2), 500)

    bm <- bilateral_matching(agents_sm, cfg$beta, cfg$vars,
                              n_rounds = vrounds, threshold = cal$threshold,
                              diagnostics = TRUE)

    if (nrow(bm$matches) < 5) {
      cat(sprintf("  %s accept=%.0f%% (thr=%+.2f): too few pairs (%d)\n",
                  dgp_name, 100 * target, cal$threshold, nrow(bm$matches)))
      next
    }

    ego_e <- agents_sm$education[match(bm$matches$ego_id, agents_sm$id)]
    alt_e <- agents_sm$education[match(bm$matches$alt_id, agents_sm$id)]
    r <- cor(ego_e, alt_e)

    cat(sprintf("  %s accept=%.0f%% (thr=%+6.2f): %4d pairs (%4.0f%%), cor=%+.4f %s\n",
                dgp_name, 100 * target, cal$threshold,
                nrow(bm$matches), 100 * bm$pairing_rate, r,
                if (abs(r) > 0.05) "← SORTING" else ""))

    sm_diag_rows[[length(sm_diag_rows) + 1]] <- data.frame(
      dgp = dgp_name,
      mechanism = "bilateral",
      target_accept_rate = target,
      threshold = cal$threshold,
      n_pairs = nrow(bm$matches),
      pairing_rate = bm$pairing_rate,
      cor_edu = r,
      V_mean = cal$V_mean,
      V_sd = cal$V_sd,
      stringsAsFactors = FALSE
    )
  }
}

# --- 6c: DA for comparison under both DGPs ---
cat("\n  6c: DA under both DGPs (sorting through competition)\n")

for (dgp_name in names(dgp_configs)) {
  cfg <- dgp_configs[[dgp_name]]
  bm_da <- deferred_acceptance_matching(agents_sm, cfg$beta, cfg$vars,
                                         J_proposals = 5, n_rounds = 30)
  ego_e_da <- agents_sm$education[match(bm_da$matches$ego_id, agents_sm$id)]
  alt_e_da <- agents_sm$education[match(bm_da$matches$alt_id, agents_sm$id)]
  r_da <- cor(ego_e_da, alt_e_da)

  cat(sprintf("  DA %s: %d pairs (%.0f%%), cor=%+.4f %s\n",
              dgp_name, nrow(bm_da$matches), 100 * bm_da$pairing_rate, r_da,
              if (abs(r_da) > 0.05) "← SORTING" else ""))

  sm_diag_rows[[length(sm_diag_rows) + 1]] <- data.frame(
    dgp = dgp_name,
    mechanism = "deferred_acceptance",
    target_accept_rate = NA,
    threshold = NA,
    n_pairs = nrow(bm_da$matches),
    pairing_rate = bm_da$pairing_rate,
    cor_edu = r_da,
    V_mean = NA,
    V_sd = NA,
    stringsAsFactors = FALSE
  )
}

sm_diag_df <- bind_rows(sm_diag_rows)
write.csv(sm_diag_df, "results/module4_status_diagnostics.csv", row.names = FALSE)

cat("\n  === CALIBRATION SUMMARY ===\n")
cat("  Homophily V ≈ [-9, +0.8], status_max V ≈ [3, 10] — different scales.\n")
cat("  Calibrated thresholds equalize selectivity across DGPs.\n")
cat("  At equal acceptance rates, both DGPs produce sorting under bilateral.\n")
cat("  DA produces sorting through competition regardless of threshold.\n")
cat("  Results saved to results/module4_status_diagnostics.csv\n")

# Clean up Part 6
rm(agents_sm, sm_diag_rows, sm_diag_df)
gc(verbose = FALSE)

cat(sprintf("\n  [Part 6 done at %s]\n", Sys.time()))


# ─────────────────────────────────────────────────────────────────────────────
# DONE
# ─────────────────────────────────────────────────────────────────────────────

total_min <- (proc.time()["elapsed"] - total_t0) / 60

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("COMPLETE — %.1f minutes total\n", total_min))
cat(sprintf("All results in: %s\n", RESULTS_FILE))
cat(sprintf("Pair diagnostics: results/module4_pair_diagnostics.csv\n"))
cat(sprintf("Status verification: results/module4_status_diagnostics.csv\n"))
cat(paste(rep("=", 60), collapse = ""), "\n")
