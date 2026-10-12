# ============================================================
# generate_appendix_tables.R
# Produces LaTeX tables for the kappa simulation appendix
# from the module CSV result files.
#
# Usage:  Rscript generate_appendix_tables.R
# Output: tables/ directory (one .tex file per table)
#
# Requires: dplyr, tidyr, readr, glue
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(glue)

# --- paths ---
results_dir <- "results"
out_dir     <- "tables"
if (!dir.exists(out_dir)) dir.create(out_dir)

fmt <- function(x, digits = 1) formatC(x, format = "f", digits = digits)
fmt_pct <- function(x) {
  s <- formatC(x, format = "f", digits = 1)
  ifelse(x > 0, paste0("$+", s, "$"), paste0("$", s, "$"))
}
fmt_est <- function(x) formatC(x, format = "f", digits = 3)


# ============================================================
# TABLE 1: Module 1 -- Group composition (oracle + bilateral)
# ============================================================
m1 <- read_csv(file.path(results_dir, "module1_results.csv"),
               show_col_types = FALSE) |>
  filter(model == "model_A") |>
  mutate(
    mech_label = case_when(
      mechanism == "oracle"          ~ "Oracle",
      grepl("bilateral", mechanism)  ~ "Bilateral"
    ),
    comp_label = case_when(
      groups == "33/33/33" ~ "Balanced (33/33/33)",
      groups == "50/30/20" ~ "Moderate (50/30/20)",
      groups == "70/20/10" ~ "Unbalanced (70/20/10)",
      groups == "85/10/5"  ~ "Extreme (85/10/5)"
    )
  ) |>
  filter(!is.na(mech_label)) |>
  select(comp_label, mech_label, parameter, rel_bias_pct, coverage_95)

m1_wide <- m1 |>
  pivot_wider(
    names_from  = parameter,
    values_from = c(rel_bias_pct, coverage_95)
  )

sink(file.path(out_dir, "tab_app_module1.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\begin{threeparttable}\n")
cat("\\caption{Module 1: Parameter recovery by group composition and matching mechanism.}\n")
cat("\\label{tab:app-module1}\n")
cat("\\begin{tabular}{ll rr rr}\n")
cat("\\toprule\n")
cat(" & & \\multicolumn{2}{c}{$\\hat{\\beta}_{\\text{edu\\_diff}}$} & \\multicolumn{2}{c}{$\\hat{\\beta}_{\\text{same\\_group}}$} \\\\\n")
cat("\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}\n")
cat("Composition & Mechanism & Bias (\\%) & Cov\\textsubscript{95} & Bias (\\%) & Cov\\textsubscript{95} \\\\\n")
cat("\\midrule\n")

comps <- c("Balanced (33/33/33)", "Moderate (50/30/20)",
           "Unbalanced (70/20/10)", "Extreme (85/10/5)")

for (i in seq_along(comps)) {
  rows <- m1_wide |> filter(comp_label == comps[i])
  for (j in seq_len(nrow(rows))) {
    r <- rows[j, ]
    comp_col <- if (j == 1) comps[i] else ""
    cat(glue("{comp_col} & {r$mech_label} & {fmt_pct(r$rel_bias_pct_edu_diff)} & {fmt(r$coverage_95_edu_diff, 2)} & {fmt_pct(r$rel_bias_pct_same_group)} & {fmt(r$coverage_95_same_group, 2)} \\\\"), "\n")
  }
  if (i < length(comps)) cat("\\addlinespace\n")
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}[flushleft]\n")
cat("\\small\n")
cat("\\item \\textit{Note.} $S = 200$. True values: $\\beta_{\\text{edu\\_diff}} = -1.5$, $\\beta_{\\text{same\\_group}} = 0.8$. Bilateral matching calibrated at 20\\% one-sided acceptance rate. Cov\\textsubscript{95} = empirical coverage of 95\\% confidence intervals.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()


# ============================================================
# TABLE 2: Module 2 -- Spatial confounding (oracle + bilateral)
# ============================================================
m2 <- read_csv(file.path(results_dir, "module2_results.csv"),
               show_col_types = FALSE)

# Select the "Both effects" experiments for edu_diff only (full + reduced)
m2_select <- m2 |>
  filter(parameter == "edu_diff",
         grepl("Both effects", experiment)) |>
  mutate(
    mech = case_when(
      grepl("bilateral", model, ignore.case = TRUE) ~ "Bilateral",
      TRUE ~ "Oracle"
    ),
    model_type = case_when(
      grepl("Reduced", model) ~ "Reduced",
      TRUE ~ "Full"
    ),
    seg_level = case_when(
      grepl("low seg", experiment)  ~ "Low",
      grepl("mod", experiment)      ~ "Moderate",
      grepl("high seg", experiment) ~ "High"
    )
  ) |>
  filter(!is.na(seg_level)) |>
  select(seg_level, mech, model_type, rel_bias_pct)

m2_wide <- m2_select |>
  pivot_wider(names_from = model_type, values_from = rel_bias_pct)

sink(file.path(out_dir, "tab_app_module2.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\begin{threeparttable}\n")
cat("\\caption{Module 2: Education homophily bias under spatial confounding.}\n")
cat("\\label{tab:app-module2}\n")
cat("\\begin{tabular}{ll rr}\n")
cat("\\toprule\n")
cat(" & & \\multicolumn{2}{c}{$\\hat{\\beta}_{\\text{edu\\_diff}}$ bias (\\%)} \\\\\n")
cat("\\cmidrule(lr){3-4}\n")
cat("Segregation & Mechanism & Full model & Reduced model \\\\\n")
cat("\\midrule\n")

segs <- c("Low", "Moderate", "High")
for (i in seq_along(segs)) {
  rows <- m2_wide |> filter(seg_level == segs[i])
  for (j in seq_len(nrow(rows))) {
    r <- rows[j, ]
    seg_col <- if (j == 1) segs[i] else ""
    cat(glue("{seg_col} & {r$mech} & {fmt_pct(r$Full)} & {fmt_pct(r$Reduced)} \\\\"), "\n")
  }
  if (i < length(segs)) cat("\\addlinespace\n")
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}[flushleft]\n")
cat("\\small\n")
cat("\\item \\textit{Note.} $S = 200$, $N = 1{,}000$. True $\\beta_{\\text{edu\\_diff}} = -1.0$, $\\beta_{\\text{log\\_distance}} = -0.5$. Full model includes edu\\_diff and log\\_distance; reduced model omits propinquity. Bilateral matching calibrated at 20\\% acceptance.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()


# ============================================================
# TABLE 3: Module 3 -- Mediation (oracle + bilateral)
# ============================================================
m3 <- read_csv(file.path(results_dir, "module3_core_stages.csv"),
               show_col_types = FALSE) |>
  filter(part == "I")

m3_edu <- m3 |> filter(parameter == "edu_diff") |>
  mutate(rel_bias_edu = bias / true_value * 100)
m3_wp  <- m3 |> filter(parameter == "same_wp") |>
  mutate(rel_bias_wp = bias / true_value * 100)

m3_tab <- m3_edu |>
  select(matching_type, stage, rel_bias_edu, mediation_pct, mediation_sd) |>
  left_join(
    m3_wp |> select(matching_type, stage, rel_bias_wp),
    by = c("matching_type", "stage")
  ) |>
  mutate(
    mech_label = case_when(
      matching_type == "one_sided"  ~ "Oracle",
      matching_type == "bilateral"  ~ "Bilateral"
    ),
    stage_label = case_when(
      stage == "A" ~ "A: No confounding",
      stage == "B" ~ "B: Confounder $\\perp$ education",
      stage == "C" ~ "C: Confounder $\\leftrightarrow$ education",
      stage == "D" ~ "D: Choice set selection"
    )
  )

sink(file.path(out_dir, "tab_app_module3.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\begin{threeparttable}\n")
cat("\\caption{Module 3: Mediation decomposition under sequential ignorability violation.}\n")
cat("\\label{tab:app-module3}\n")
cat("\\begin{tabular}{ll rr rr}\n")
cat("\\toprule\n")
cat(" & & \\multicolumn{2}{c}{Parameter bias (\\%)} & & \\\\\n")
cat("\\cmidrule(lr){3-4}\n")
cat("Stage & Mechanism & $\\hat{\\beta}_{\\text{edu}}$ & $\\hat{\\beta}_{\\text{wp}}$ & Med.\\ \\% & Med.\\ SD \\\\\n")
cat("\\midrule\n")

stages <- c("A", "B", "C", "D")
stage_labels <- c("A: No confounding",
                  "B: Confounder $\\perp$ education",
                  "C: Confounder $\\leftrightarrow$ education",
                  "D: Choice set selection")

for (i in seq_along(stages)) {
  rows <- m3_tab |> filter(stage == stages[i])
  for (j in seq_len(nrow(rows))) {
    r <- rows[j, ]
    st_col <- if (j == 1) stage_labels[i] else ""
    cat(glue("{st_col} & {r$mech_label} & {fmt_pct(r$rel_bias_edu)} & {fmt_pct(r$rel_bias_wp)} & {fmt(r$mediation_pct)} & {fmt(r$mediation_sd)} \\\\"), "\n")
  }
  if (i < length(stages)) cat("\\addlinespace\n")
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}[flushleft]\n")
cat("\\small\n")
cat("\\item \\textit{Note.} $S = 200$. True values: $\\beta_{\\text{edu}} = -1.0$, $\\beta_{\\text{wp}} = 1.5$. Stage~A: no confounder. Stage~B: ``openness'' confounds workplace sorting and partner choice, uncorrelated with education. Stage~C: openness correlated with education ($\\rho = 0.3$). Stage~D: choice set selection (oracle only). Med.\\ \\% = share of the education--partner association mediated through the workplace. Bilateral calibrated at 20\\% acceptance.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()


# ============================================================
# TABLE 4: Module 4 -- Mechanism comparison (oracle/bilateral/DA)
# ============================================================
m4 <- read_csv(file.path(results_dir, "module4_results.csv"),
               show_col_types = FALSE) |>
  filter(part == "P1_core")

# Build the table for homophily and status max DGPs
m4_h <- m4 |> filter(dgp == "homophily")
m4_s <- m4 |> filter(dgp == "status_max")
m4_m <- m4 |> filter(dgp == "mixed")

make_m4_row <- function(df, dgp_label, mech_filter, edu_param = "edu_diff") {
  row <- df |> filter(mechanism == mech_filter)
  edu <- row |> filter(parameter == edu_param)
  grp <- row |> filter(parameter == "same_group")
  mech_label <- case_when(
    mech_filter == "one_sided_oracle"   ~ "Oracle",
    mech_filter == "bilateral"          ~ "Bilateral",
    mech_filter == "deferred_acceptance" ~ "DA"
  )
  glue(" & {mech_label} & ${fmt_est(edu$mean_est)}$ & {fmt_pct(edu$rel_bias_pct)} & ${fmt_est(grp$mean_est)}$ & {fmt_pct(grp$rel_bias_pct)} & {fmt(edu$mean_pair_rate, 2)} \\\\")
}

sink(file.path(out_dir, "tab_app_module4.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\begin{threeparttable}\n")
cat("\\caption{Module 4: Parameter recovery by matching mechanism and DGP.}\n")
cat("\\label{tab:app-module4}\n")
cat("\\begin{tabular}{ll rr rr r}\n")
cat("\\toprule\n")
cat(" & & \\multicolumn{2}{c}{$\\hat{\\beta}_{\\text{edu}}$} & \\multicolumn{2}{c}{$\\hat{\\beta}_{\\text{group}}$} & \\\\\n")
cat("\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}\n")
cat("DGP & Mechanism & $\\bar{\\hat{\\beta}}$ & Bias (\\%) & $\\bar{\\hat{\\beta}}$ & Bias (\\%) & Pair rate \\\\\n")
cat("\\midrule\n")

cat("\\multirow{3}{*}{Homophily}", make_m4_row(m4_h, "Homophily", "one_sided_oracle"), "\n")
cat(make_m4_row(m4_h, "Homophily", "bilateral"), "\n")
cat(make_m4_row(m4_h, "Homophily", "deferred_acceptance"), "\n")
cat("\\addlinespace\n")

cat("\\multirow{3}{*}{Status max.}", make_m4_row(m4_s, "Status max.", "one_sided_oracle", "alt_edu"), "\n")
cat(make_m4_row(m4_s, "Status max.", "bilateral", "alt_edu"), "\n")
cat(make_m4_row(m4_s, "Status max.", "deferred_acceptance", "alt_edu"), "\n")

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}[flushleft]\n")
cat("\\small\n")
cat("\\item \\textit{Note.} $S = 200$. Oracle = one-sided choice (the process conditional logit assumes). Bilateral = random encounters with calibrated threshold (20\\% acceptance). DA = deferred acceptance. Pair rate = mean partners per agent.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()


cat("Tables written to", out_dir, "\n")
