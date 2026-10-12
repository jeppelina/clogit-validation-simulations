# =============================================================================
# estimate.R -- choice sets and conditional logit, as in the papers
# =============================================================================
# For each partnered ego: the actual partner plus J - 1 alternatives drawn
# uniformly from the rest of the opposite-sex pool. One model per ego sex.
# =============================================================================

library(survival)

build_choice_sets <- function(ego_ids, partner_ids, ego_df, alt_df, J = 100) {
  n <- length(ego_ids); N_alt <- nrow(alt_df)
  alts <- t(vapply(seq_len(n), function(i) {
    c(partner_ids[i], sample(setdiff(seq_len(N_alt), partner_ids[i]), J - 1))
  }, integer(J)))
  ego <- rep(ego_ids, each = J)
  alt <- as.vector(t(alts))
  data.frame(
    chooser = ego,
    chosen = rep(c(1L, rep(0L, J - 1)), n),
    same_group = as.integer(ego_df$group[ego] == alt_df$group[alt]),
    edu_diff = abs(ego_df$edu[ego] - alt_df$edu[alt])
  )
}

fit_clogit <- function(cs, terms = c("same_group", "edu_diff")) {
  f <- as.formula(paste("chosen ~", paste(terms, collapse = " + "),
                        "+ strata(chooser)"))
  fit <- clogit(f, data = cs, method = "breslow")
  data.frame(term = terms, est = unname(coef(fit)[terms]),
             se = unname(sqrt(diag(vcov(fit)))[match(terms, names(coef(fit)))]))
}
