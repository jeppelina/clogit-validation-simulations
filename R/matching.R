# =============================================================================
# matching.R: how partners are formed
# =============================================================================
# one_sided():  each person picks the best option in the whole other-sex pool
#               (or stays single). The chosen person has no say and can be
#               chosen by many. This is the process conditional logit assumes.
# gale_shapley(): both partners must prefer the union to staying single.
#               Men propose in order of preference; each woman keeps the best
#               offer so far. The result is a stable matching: no man and
#               woman would both rather be with each other.
# =============================================================================

#' One-sided choice. U: [chooser, alternative] utilities; out: value of
#' staying single. Returns the chosen alternative per chooser (NA = single).
one_sided <- function(U, out) {
  best <- max.col(U, ties.method = "first")
  ifelse(U[cbind(seq_len(nrow(U)), best)] > out, best, NA_integer_)
}

#' Men-proposing deferred acceptance (Gale-Shapley) with an outside option.
#' Returns a data.frame of couples (man, woman).
gale_shapley <- function(ut) {
  n_m <- nrow(ut$U_m); n_w <- ncol(ut$U_m)
  lists <- lapply(seq_len(n_m), function(m) {
    acceptable <- which(ut$U_m[m, ] > ut$out_m[m])
    acceptable[order(ut$U_m[m, acceptable], decreasing = TRUE)]
  })
  next_i <- rep(1L, n_m)
  held <- rep(NA_integer_, n_w)                    # woman -> man she holds
  queue <- which(lengths(lists) > 0)
  while (length(queue) > 0) {
    m <- queue[1]; queue <- queue[-1]
    if (next_i[m] > length(lists[[m]])) next       # no one left: stays single
    w <- lists[[m]][next_i[m]]; next_i[m] <- next_i[m] + 1L
    u <- ut$U_w[m, w]
    if (u <= ut$out_w[w]) { queue <- c(queue, m); next }
    h <- held[w]
    if (is.na(h)) held[w] <- m
    else if (u > ut$U_w[h, w]) { held[w] <- m; queue <- c(queue, h) }
    else queue <- c(queue, m)
  }
  ok <- !is.na(held)
  data.frame(man = held[ok], woman = which(ok))
}

#' Number of blocking pairs: a man and a woman who both strictly prefer each
#' other to their outcome. Zero for a stable matching.
count_blocking_pairs <- function(ut, couples) {
  cur_m <- ut$out_m; cur_w <- ut$out_w
  cur_m[couples$man] <- ut$U_m[cbind(couples$man, couples$woman)]
  cur_w[couples$woman] <- ut$U_w[cbind(couples$man, couples$woman)]
  block <- (ut$U_m > cur_m) & sweep(ut$U_w, 2, cur_w, ">")
  block[cbind(couples$man, couples$woman)] <- FALSE
  sum(block)
}

#' Find the value of staying single (u0) that makes a target share of people
#' pair up under Gale-Shapley. make_ut(u0) must build a market from fixed
#' random numbers; bisection over u0.
calibrate_u0 <- function(make_ut, target = 0.85, lower = -5, upper = 15, iters = 14) {
  rate_at <- function(u0) mean(vapply(1:2, function(k) {
    set.seed(9000 + k); ut <- make_ut(u0)
    nrow(gale_shapley(ut)) / nrow(ut$U_m)
  }, numeric(1)))
  for (i in seq_len(iters)) {
    mid <- (lower + upper) / 2
    if (rate_at(mid) > target) lower <- mid else upper <- mid
  }
  (lower + upper) / 2
}
