# =============================================================================
# population.R -- stylised two-sex partner market
# =============================================================================
# One market: N men and N women. Each person has an ancestry group and a
# continuous education score. Utilities are held as N x N matrices indexed
# [man, woman], one per side, with taste shocks drawn ONCE and kept fixed for
# the whole matching process.
# =============================================================================

# seg: residential segregation, 0 to 1, or NULL for no space at all. With
# space, everyone lives on the unit square: with probability seg near their
# group's centre (sd 0.08), otherwise anywhere at random.
GROUP_CENTRES <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.5, 0.75),
                       c(0.25, 0.75), c(0.75, 0.75))

# shares_w / edu_shift_w let the women's side differ from the men's (Module
# B). By default both sexes are drawn the same way.
# layout / s (Module A v3, two groups): residential layouts on the unit square.
#   "mixed"    everyone uniform
#   "sides"    with probability s a person lives on their group's own side:
#              the majority on x in [0, share_majority], the minority on the
#              rest, so both groups have the same density when fully separated
#   "enclaves" the minority, with probability s, lives in one of 3 small
#              clusters (sd 0.05); the majority, with probability s, avoids
#              the clusters (radius 0.15). Same centres for both sexes.
#   Otherwise (probability 1 - s) a person lives anywhere.
make_population <- function(N, shares = c(0.7, 0.2, 0.1), edu_sd = 1,
                            seg = NULL, shares_w = shares, edu_shift_w = 0,
                            layout = NULL, s = 0) {
  draw_side <- function(sh, shift = 0) {
    sizes <- round(N * sh)
    sizes[1] <- N - sum(sizes[-1])
    d <- data.frame(id = seq_len(N),
                    group = rep(seq_along(sh), times = sizes),
                    edu = rnorm(N, shift, edu_sd))
    if (!is.null(seg)) {
      near <- runif(N) < seg
      cx <- GROUP_CENTRES[d$group, 1]; cy <- GROUP_CENTRES[d$group, 2]
      d$x <- ifelse(near, pmin(pmax(rnorm(N, cx, 0.08), 0), 1), runif(N))
      d$y <- ifelse(near, pmin(pmax(rnorm(N, cy, 0.08), 0), 1), runif(N))
    }
    d
  }
  pop <- list(men = draw_side(shares), women = draw_side(shares_w, edu_shift_w))
  if (!is.null(layout)) pop <- place_two_groups(pop, layout, s, shares[1])
  pop
}

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
      ma <- which(own & d$group == 1)
      near <- function(i) any(sqrt((x[i] - ENCLAVES[, 1])^2 + (y[i] - ENCLAVES[, 2])^2) < 0.15)
      for (i in ma) for (t in 1:50) { if (!near(i)) break; x[i] <- runif(1); y[i] <- runif(1) }
    }
    d$x <- x; d$y <- y
    d
  }
  list(men = place(pop$men), women = place(pop$women))
}

# Institutions (Module C). K institutions, each with an education level c_k
# spread over the education range (and a trait level for the unobserved
# trait t, if present). People join institution k with probability
# proportional to exp(-sort * distance to (c_k, t_k)), so institutions sort on
# education (and t). trait = TRUE adds t ~ N(0, 1) to everyone.
add_institutions <- function(pop, K = 20, sort = 2, trait = FALSE) {
  c_edu <- qnorm(seq(0.5 / K, 1 - 0.5 / K, length.out = K))
  c_t <- if (trait) sample(c_edu) else rep(0, K)
  assign_side <- function(d) {
    if (trait) d$t <- rnorm(nrow(d))
    tt <- if (trait) d$t else rep(0, nrow(d))
    dist <- abs(outer(d$edu, c_edu, "-")) + abs(outer(tt, c_t, "-"))
    w <- exp(-sort * dist)
    d$inst <- apply(w, 1, function(p) sample.int(K, 1, prob = p))
    d
  }
  list(men = assign_side(pop$men), women = assign_side(pop$women))
}

ln_dist <- function(x1, y1, x2, y2) log(sqrt((x1 - x2)^2 + (y1 - y2)^2) + 0.01)

# Systematic utility for every man-woman pair, from one side's coefficients.
# beta: named vector with any of same_group, edu_diff, ln_dist, alt_edu.
# alt_edu is the partner's own education (status), so it depends on whose
# utility this is: side = "men" (man values the woman) or "women".
pair_V <- function(pop, beta, side = "men") {
  V <- matrix(0, nrow(pop$men), nrow(pop$women))
  if ("same_group" %in% names(beta))
    V <- V + beta[["same_group"]] * outer(pop$men$group, pop$women$group, "==")
  if ("edu_diff" %in% names(beta))
    V <- V + beta[["edu_diff"]] * abs(outer(pop$men$edu, pop$women$edu, "-"))
  if ("ln_dist" %in% names(beta))
    V <- V + beta[["ln_dist"]] * log(sqrt(outer(pop$men$x, pop$women$x, "-")^2 +
                                          outer(pop$men$y, pop$women$y, "-")^2) + 0.01)
  if ("same_inst" %in% names(beta))
    V <- V + beta[["same_inst"]] * outer(pop$men$inst, pop$women$inst, "==")
  if ("t_diff" %in% names(beta))
    V <- V + beta[["t_diff"]] * abs(outer(pop$men$t, pop$women$t, "-"))
  if ("alt_edu" %in% names(beta))
    V <- V + beta[["alt_edu"]] * (if (side == "men")
      matrix(pop$women$edu, nrow(pop$men), nrow(pop$women), byrow = TRUE) else
      matrix(pop$men$edu, nrow(pop$men), nrow(pop$women)))
  V
}

rgumbel_matrix <- function(n, m) matrix(-log(-log(runif(n * m))), n, m)

# Full utilities with fixed shocks. Both matrices are indexed [man, woman]:
#   U_m[m, w] = what man m gets from woman w
#   U_w[m, w] = what woman w gets from man m
# Outside options: u0 + Gumbel, one fixed draw per person. u0 = -Inf means
# nobody prefers to stay single.
make_utilities <- function(pop, beta_m, beta_w, u0 = -Inf) {
  N_m <- nrow(pop$men); N_w <- nrow(pop$women)
  list(
    U_m = pair_V(pop, beta_m, "men") + rgumbel_matrix(N_m, N_w),
    U_w = pair_V(pop, beta_w, "women") + rgumbel_matrix(N_m, N_w),
    out_m = u0 + -log(-log(runif(N_m))),
    out_w = u0 + -log(-log(runif(N_w)))
  )
}
