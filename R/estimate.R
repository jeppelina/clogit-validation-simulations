# =============================================================================
# estimate.R -- choice sets and conditional logit, as in the papers
# =============================================================================
# For each partnered ego: the actual partner plus J - 1 alternatives from the
# rest of the opposite-sex pool. One model per ego sex.
#   build_choice_sets():       uniform sampling.
#   build_choice_sets_strat(): stratified sampling with a fixed number of
#     same-group alternatives (oversampled), plus the McFadden-Manski offset
#     log(N_s / n_s) for each alternative's stratum s, where N_s is the
#     stratum's pool size and n_s its count in the choice set.
# =============================================================================

library(survival)

pair_vars <- function(ego, alt, ego_df, alt_df) {
  d <- data.frame(
    same_group = as.integer(ego_df$group[ego] == alt_df$group[alt]),
    edu_diff = abs(ego_df$edu[ego] - alt_df$edu[alt]),
    alt_edu = alt_df$edu[alt])
  # coarse education: same tercile of N(0, 1) education (a stand-in for
  # categorical homogamy measures such as Paper 4's Edudiff3)
  terc <- function(e) findInterval(e, qnorm(c(1, 2) / 3))
  d$same_edu3 <- as.integer(terc(ego_df$edu[ego]) == terc(alt_df$edu[alt]))
  # alter-group indicators (reference: group 1), as in Paper 3's alter-category
  # specification
  for (k in setdiff(sort(unique(alt_df$group)), 1))
    d[[paste0("alt_g", k)]] <- as.integer(alt_df$group[alt] == k)
  if (!is.null(ego_df$inst))
    d$same_inst <- as.integer(ego_df$inst[ego] == alt_df$inst[alt])
  if (!is.null(ego_df$t))
    d$t_diff <- abs(ego_df$t[ego] - alt_df$t[alt])
  if (!is.null(ego_df$x))
    d$ln_dist <- ln_dist(ego_df$x[ego], ego_df$y[ego], alt_df$x[alt], alt_df$y[alt])
  d
}

build_choice_sets <- function(ego_ids, partner_ids, ego_df, alt_df, J = 100) {
  n <- length(ego_ids); N_alt <- nrow(alt_df)
  alts <- t(vapply(seq_len(n), function(i) {
    c(partner_ids[i], sample(setdiff(seq_len(N_alt), partner_ids[i]), J - 1))
  }, integer(J)))
  ego <- rep(ego_ids, each = J)
  alt <- as.vector(t(alts))
  cbind(data.frame(chooser = ego, chosen = rep(c(1L, rep(0L, J - 1)), n)),
        pair_vars(ego, alt, ego_df, alt_df))
}

# n_same of the J alternatives come from the ego's own group (counting the
# partner if same-group), the rest from other groups.
build_choice_sets_strat <- function(ego_ids, partner_ids, ego_df, alt_df,
                                    J = 100, n_same = 50) {
  N_alt <- nrow(alt_df)
  rows <- lapply(seq_along(ego_ids), function(i) {
    e <- ego_ids[i]; p <- partner_ids[i]; g <- ego_df$group[e]
    same_pool <- which(alt_df$group == g); other_pool <- which(alt_df$group != g)
    p_same <- alt_df$group[p] == g
    k_same <- n_same - p_same; k_other <- (J - n_same) - !p_same
    s_draw <- setdiff(same_pool, p); o_draw <- setdiff(other_pool, p)
    alts <- c(p, s_draw[sample.int(length(s_draw), min(k_same, length(s_draw)))],
              o_draw[sample.int(length(o_draw), min(k_other, length(o_draw)))])
    same <- alt_df$group[alts] == g
    n_s <- sum(same); n_o <- sum(!same)
    data.frame(chooser = e, alt = alts, chosen = c(1L, rep(0L, length(alts) - 1)),
               offset_mm = ifelse(same, log(length(same_pool) / n_s),
                                  log(length(other_pool) / n_o)))
  })
  d <- do.call(rbind, rows)
  cbind(d[, c("chooser", "chosen", "offset_mm")],
        pair_vars(d$chooser, d$alt, ego_df, alt_df))
}

fit_clogit <- function(cs, terms = c("same_group", "edu_diff"), offset = FALSE) {
  rhs <- paste(terms, collapse = " + ")
  if (offset) rhs <- paste(rhs, "+ offset(offset_mm)")
  f <- as.formula(paste("chosen ~", rhs, "+ strata(chooser)"))
  fit <- clogit(f, data = cs, method = "breslow")
  data.frame(term = terms, est = unname(coef(fit)[terms]),
             se = unname(sqrt(diag(vcov(fit)))[match(terms, names(coef(fit)))]))
}

# KHB decomposition (Karlson, Holm and Breen 2012) for one mediator.
# Full model: focal + mediator + controls. Reduced model: focal + residualised
# mediator + controls, where the mediator is residualised on focal and
# controls by OLS over all choice-set rows. Both models are on the same scale,
# so the change in the focal coefficient is the mediated part.
khb_share <- function(cs, focal, mediator, controls = character(0)) {
  rhs_x <- paste(c(focal, controls), collapse = " + ")
  cs$.z_res <- resid(lm(as.formula(paste(mediator, "~", rhs_x)), data = cs))
  f_full <- clogit(as.formula(paste("chosen ~", rhs_x, "+", mediator, "+ strata(chooser)")),
                   data = cs, method = "breslow")
  f_red <- clogit(as.formula(paste("chosen ~", rhs_x, "+ .z_res + strata(chooser)")),
                  data = cs, method = "breslow")
  b_full <- coef(f_full)[[focal]]; b_red <- coef(f_red)[[focal]]
  data.frame(b_reduced = b_red, b_full = b_full, b_mediator = coef(f_full)[[mediator]],
             share = (b_red - b_full) / b_red,
             cor_med_focal = cor(cs[[mediator]], cs[[focal]]))
}

# Predicted probability mass on alternatives with flag = 1 (e.g. same group),
# averaged over choosers, at coefficient vector b (named by term).
predicted_share <- function(cs, b, flag = "same_group") {
  eta <- as.matrix(cs[, names(b)]) %*% b
  eta <- eta - ave(eta, cs$chooser, FUN = max)
  p <- exp(eta); p <- p / ave(p, cs$chooser, FUN = sum)
  mean(tapply(p * cs[[flag]], cs$chooser, sum))
}
