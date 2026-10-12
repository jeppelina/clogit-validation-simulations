# =============================================================================
# population.R -- stylised two-sex partner market
# =============================================================================
# One market: N men and N women. Each person has an ancestry group and a
# continuous education score. Utilities are held as N x N matrices indexed
# [man, woman], one per side, with taste shocks drawn ONCE and kept fixed for
# the whole matching process.
# =============================================================================

make_population <- function(N, shares = c(0.7, 0.2, 0.1), edu_sd = 1) {
  draw_side <- function() {
    sizes <- round(N * shares)
    sizes[1] <- N - sum(sizes[-1])
    data.frame(id = seq_len(N),
               group = rep(seq_along(shares), times = sizes),
               edu = rnorm(N, 0, edu_sd))
  }
  list(men = draw_side(), women = draw_side())
}

# Systematic utility for every man-woman pair, from one side's coefficients.
# beta: named vector with any of same_group, edu_diff.
pair_V <- function(pop, beta) {
  V <- matrix(0, nrow(pop$men), nrow(pop$women))
  if (!is.null(beta[["same_group"]]))
    V <- V + beta[["same_group"]] * outer(pop$men$group, pop$women$group, "==")
  if (!is.null(beta[["edu_diff"]]))
    V <- V + beta[["edu_diff"]] * abs(outer(pop$men$edu, pop$women$edu, "-"))
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
    U_m = pair_V(pop, beta_m) + rgumbel_matrix(N_m, N_w),
    U_w = pair_V(pop, beta_w) + rgumbel_matrix(N_m, N_w),
    out_m = u0 + -log(-log(runif(N_m))),
    out_w = u0 + -log(-log(runif(N_w)))
  )
}
