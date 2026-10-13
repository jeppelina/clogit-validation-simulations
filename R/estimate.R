# =============================================================================
# estimate.R: choice sets, conditional logit, KHB, predicted shares
# =============================================================================
# Each partnered ego's choice set is the actual partner plus J - 1
# alternatives drawn at random from the rest of the other sex.
# =============================================================================

library(survival)

#' Pair variables for ego rows `ego` (in ego_df) and alternatives `alt`
#' (in alt_df).
pair_vars <- function(ego, alt, ego_df, alt_df) {
  d <- data.frame(same_group = as.integer(ego_df$group[ego] == alt_df$group[alt]),
                  edu_diff = abs(ego_df$edu[ego] - alt_df$edu[alt]),
                  alt_edu = alt_df$edu[alt],
                  alt_minority = as.integer(alt_df$group[alt] == 2))
  if (!is.null(ego_df$x)) d$ln_dist <- ln_dist(ego_df$x[ego], ego_df$y[ego], alt_df$x[alt], alt_df$y[alt])
  if (!is.null(ego_df$setting)) d$same_setting <- as.integer(ego_df$setting[ego] == alt_df$setting[alt])
  if (!is.null(ego_df$t)) d$t_diff <- abs(ego_df$t[ego] - alt_df$t[alt])
  d
}

build_choice_sets <- function(ego_ids, partner_ids, ego_df, alt_df, J = 100) {
  n <- length(ego_ids); N_alt <- nrow(alt_df)
  alts <- t(vapply(seq_len(n), function(i)
    c(partner_ids[i], sample(setdiff(seq_len(N_alt), partner_ids[i]), J - 1)), integer(J)))
  ego <- rep(ego_ids, each = J); alt <- as.vector(t(alts))
  cbind(data.frame(chooser = ego, chosen = rep(c(1L, rep(0L, J - 1)), n)),
        pair_vars(ego, alt, ego_df, alt_df))
}

#' Conditional logit; returns term, estimate, standard error.
fit_clogit <- function(cs, terms) {
  f <- as.formula(paste("chosen ~", paste(terms, collapse = " + "), "+ strata(chooser)"))
  fit <- clogit(f, data = cs, method = "breslow")
  data.frame(term = terms, est = unname(coef(fit)[terms]),
             se = unname(sqrt(diag(vcov(fit)))[match(terms, names(coef(fit)))]))
}

#' KHB mediation (Karlson, Holm and Breen 2012). Full model: focal + mediator
#' + controls. Reduced model: focal + mediator residualised on focal and
#' controls + controls, which keeps both models on the same scale.
#' share = (b_reduced - b_full) / b_reduced.
khb_share <- function(cs, focal, mediator, controls = character(0)) {
  rhs <- paste(c(focal, controls), collapse = " + ")
  cs$.m_res <- resid(lm(as.formula(paste(mediator, "~", rhs)), data = cs))
  full <- clogit(as.formula(paste("chosen ~", rhs, "+", mediator, "+ strata(chooser)")), data = cs, method = "breslow")
  red <- clogit(as.formula(paste("chosen ~", rhs, "+ .m_res + strata(chooser)")), data = cs, method = "breslow")
  b_full <- coef(full)[[focal]]; b_red <- coef(red)[[focal]]
  data.frame(b_reduced = b_red, b_full = b_full, b_mediator = coef(full)[[mediator]],
             share = (b_red - b_full) / b_red)
}

#' Predicted share of choices on alternatives with flag = 1 (e.g. same
#' group), averaged over choosers, at coefficients b. At the estimates this
#' equals the observed share (the score equation for that term).
predicted_share <- function(cs, b, flag = "same_group") {
  eta <- as.matrix(cs[, names(b)]) %*% b
  eta <- eta - ave(eta, cs$chooser, FUN = max)
  p <- exp(eta); p <- p / ave(p, cs$chooser, FUN = sum)
  mean(tapply(p * cs[[flag]], cs$chooser, sum))
}
