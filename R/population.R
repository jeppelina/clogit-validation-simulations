# =============================================================================
# population.R: a stylised two-sex partner market
# =============================================================================
# A market has N men and N women. Each person has a group and a continuous
# education score, and optionally a location, a setting (e.g. a workplace)
# and an unobserved trait. Utilities are N x N matrices indexed
# [man, woman], one for each side, with random tastes drawn once.
# =============================================================================

#' Draw a population.
#' @param N people per sex
#' @param shares group shares for men (group 1 = majority)
#' @param shares_w group shares for women (default: as for men)
#' @param layout NULL (no locations), "mixed", "sides" or "enclaves"
#' @param s segregation: probability that a person lives in their group's area
make_population <- function(N, shares = c(0.8, 0.2), shares_w = shares,
                            layout = NULL, s = 0) {
  draw_side <- function(sh) {
    sizes <- round(N * sh)
    sizes[1] <- N - sum(sizes[-1])
    data.frame(id = seq_len(N), group = rep(seq_along(sh), times = sizes),
               edu = rnorm(N))
  }
  pop <- list(men = draw_side(shares), women = draw_side(shares_w))
  if (!is.null(layout)) pop <- place_two_groups(pop, layout, s, shares[1])
  pop
}

# Residential layouts on the unit square (two groups).
#   mixed:    everyone lives anywhere.
#   sides:    with probability s a person lives on their group's side; the
#             majority's side is [0, share_majority] in x, so both groups
#             have the same density when fully separated.
#   enclaves: with probability s a minority person lives in one of three
#             small clusters (sd 0.05) and a majority person avoids them
#             (radius 0.15).
ENCLAVES <- rbind(c(0.25, 0.70), c(0.70, 0.75), c(0.55, 0.25))

place_two_groups <- function(pop, layout, s, share_majority) {
  place <- function(d) {
    n <- nrow(d); x <- runif(n); y <- runif(n)
    own <- runif(n) < s
    if (layout == "sides") {
      maj <- d$group == 1
      x[own & maj] <- runif(sum(own & maj), 0, share_majority)
      x[own & !maj] <- runif(sum(own & !maj), share_majority, 1)
    } else if (layout == "enclaves") {
      mi <- which(own & d$group == 2)
      k <- sample.int(nrow(ENCLAVES), length(mi), replace = TRUE)
      x[mi] <- pmin(pmax(rnorm(length(mi), ENCLAVES[k, 1], 0.05), 0), 1)
      y[mi] <- pmin(pmax(rnorm(length(mi), ENCLAVES[k, 2], 0.05), 0), 1)
      near <- function(i) any(sqrt((x[i] - ENCLAVES[, 1])^2 + (y[i] - ENCLAVES[, 2])^2) < 0.15)
      for (i in which(own & d$group == 1))
        for (t in 1:50) { if (!near(i)) break; x[i] <- runif(1); y[i] <- runif(1) }
    }
    d$x <- x; d$y <- y
    d
  }
  list(men = place(pop$men), women = place(pop$women))
}

#' Add settings (e.g. workplaces) and an unobserved trait t ~ N(0, 1).
#' K settings each have an education level c_k and a trait level d_k. A
#' person joins setting k with probability proportional to
#' exp(-sort * (|edu - c_k| + |t - d_k|)), so settings sort on both.
add_settings <- function(pop, K = 10, sort = 5) {
  c_edu <- qnorm(seq(0.5 / K, 1 - 0.5 / K, length.out = K))
  d_t <- sample(c_edu)
  assign_side <- function(d) {
    d$t <- rnorm(nrow(d))
    w <- exp(-sort * (abs(outer(d$edu, c_edu, "-")) + abs(outer(d$t, d_t, "-"))))
    d$setting <- apply(w, 1, function(p) sample.int(K, 1, prob = p))
    d
  }
  list(men = assign_side(pop$men), women = assign_side(pop$women))
}

ln_dist <- function(x1, y1, x2, y2) log(sqrt((x1 - x2)^2 + (y1 - y2)^2) + 0.01)

#' Systematic utility for every man-woman pair from one side's tastes.
#' beta: named vector with any of same_group, edu_diff (|edu difference|),
#' ln_dist, alt_edu (the partner's education), same_setting, t_diff.
#' side: whose utility it is ("men": a man valuing a woman, or "women").
pair_V <- function(pop, beta, side = "men") {
  m <- pop$men; w <- pop$women
  has <- function(term) term %in% names(beta)
  V <- matrix(0, nrow(m), nrow(w))
  if (has("same_group")) V <- V + beta[["same_group"]] * outer(m$group, w$group, "==")
  if (has("edu_diff")) V <- V + beta[["edu_diff"]] * abs(outer(m$edu, w$edu, "-"))
  if (has("ln_dist")) V <- V + beta[["ln_dist"]] * log(sqrt(outer(m$x, w$x, "-")^2 + outer(m$y, w$y, "-")^2) + 0.01)
  if (has("same_setting")) V <- V + beta[["same_setting"]] * outer(m$setting, w$setting, "==")
  if (has("t_diff")) V <- V + beta[["t_diff"]] * abs(outer(m$t, w$t, "-"))
  if (has("alt_edu")) V <- V + beta[["alt_edu"]] * (if (side == "men")
    matrix(w$edu, nrow(m), nrow(w), byrow = TRUE) else matrix(m$edu, nrow(m), nrow(w)))
  V
}

rgumbel_matrix <- function(n, m) matrix(-log(-log(runif(n * m))), n, m)

#' Utilities with random tastes drawn once. Both matrices are [man, woman]:
#'   U_m[m, w] = what man m gets from woman w
#'   U_w[m, w] = what woman w gets from man m
#' Value of staying single: u0 + Gumbel, one draw per person (u0 = -Inf:
#' nobody prefers to stay single).
make_utilities <- function(pop, beta_m, beta_w, u0 = -Inf) {
  n_m <- nrow(pop$men); n_w <- nrow(pop$women)
  list(U_m = pair_V(pop, beta_m, "men") + rgumbel_matrix(n_m, n_w),
       U_w = pair_V(pop, beta_w, "women") + rgumbel_matrix(n_m, n_w),
       out_m = u0 - log(-log(runif(n_m))),
       out_w = u0 - log(-log(runif(n_w))))
}
