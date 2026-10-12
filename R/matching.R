# =============================================================================
# matching.R -- how partners are formed
# =============================================================================
# Two mechanisms.
#
# oracle_choice(): each ego picks the best option in the whole opposite-sex
#   pool, or stays single. No capacity limits, the other side has no say.
#   This is exactly the process conditional logit assumes. Run per sex.
#
# deferred_acceptance(): people meet k random others (a fixed meeting
#   network), have fixed tastes, can stay single, and among those who met,
#   proposers and receivers compete until no proposer has anyone left to try.
#   k = everyone gives textbook Gale-Shapley. The result is stable with
#   respect to the meeting network: no two people who MET would both rather
#   be with each other than with their outcome.
# =============================================================================

oracle_choice <- function(U_ego, out_ego) {
  # U_ego: [ego, alter] utilities. Returns the chosen alter per ego, NA = single.
  best <- max.col(U_ego, ties.method = "first")
  best_u <- U_ego[cbind(seq_len(nrow(U_ego)), best)]
  ifelse(best_u > out_ego, best, NA_integer_)
}

# Random bipartite meeting network: each man meets k women (women's number
# of meetings then varies around k). Meetings are mutual. k >= N: everyone.
meeting_network <- function(N_m, N_w, k) {
  if (k >= N_w) return(matrix(TRUE, N_m, N_w))
  M <- matrix(FALSE, N_m, N_w)
  for (m in seq_len(N_m)) M[m, sample.int(N_w, k)] <- TRUE
  M
}

# Deferred acceptance on a meeting network.
#   U_prop[p, r]: proposer p's utility of receiver r
#   U_recv[p, r]: receiver r's utility of proposer p
# Returns partner_of_proposer (receiver index or NA).
deferred_acceptance <- function(U_prop, U_recv, out_prop, out_recv, met) {
  n_p <- nrow(U_prop); n_r <- ncol(U_prop)
  # Each proposer's list: met receivers they prefer to staying single, best first.
  lists <- lapply(seq_len(n_p), function(p) {
    cand <- which(met[p, ] & U_prop[p, ] > out_prop[p])
    cand[order(U_prop[p, cand], decreasing = TRUE)]
  })
  next_i <- rep(1L, n_p)
  held <- rep(NA_integer_, n_r)      # receiver -> proposer currently held
  queue <- which(lengths(lists) > 0)
  while (length(queue) > 0) {
    p <- queue[1]; queue <- queue[-1]
    if (next_i[p] > length(lists[[p]])) next          # list exhausted: single
    r <- lists[[p]][next_i[p]]; next_i[p] <- next_i[p] + 1L
    u <- U_recv[p, r]
    if (u <= out_recv[r]) { queue <- c(queue, p); next }   # r rejects outright
    h <- held[r]
    if (is.na(h)) {
      held[r] <- p
    } else if (u > U_recv[h, r]) {
      held[r] <- p; queue <- c(queue, h)                   # h is displaced
    } else {
      queue <- c(queue, p)
    }
  }
  partner <- rep(NA_integer_, n_p)
  ok <- !is.na(held)
  partner[held[ok]] <- which(ok)
  partner
}

# Wrapper: returns a data.frame of couples (man, woman) under either side
# proposing.
match_market <- function(ut, k, proposer = c("men", "women")) {
  proposer <- match.arg(proposer)
  N_m <- nrow(ut$U_m); N_w <- ncol(ut$U_m)
  met <- meeting_network(N_m, N_w, k)
  if (proposer == "men") {
    pm <- deferred_acceptance(ut$U_m, ut$U_w, ut$out_m, ut$out_w, met)
    keep <- !is.na(pm)
    data.frame(man = which(keep), woman = pm[keep])
  } else {
    pw <- deferred_acceptance(t(ut$U_w), t(ut$U_m), ut$out_w, ut$out_m, t(met))
    keep <- !is.na(pw)
    data.frame(man = pw[keep], woman = which(keep))
  }
}

# Stability check against the meeting network: returns the number of
# blocking pairs (met, both strictly prefer each other to their outcome).
count_blocking_pairs <- function(ut, couples, met) {
  N_m <- nrow(ut$U_m); N_w <- ncol(ut$U_m)
  cur_m <- ut$out_m; cur_w <- ut$out_w
  cur_m[couples$man] <- ut$U_m[cbind(couples$man, couples$woman)]
  cur_w[couples$woman] <- ut$U_w[cbind(couples$man, couples$woman)]
  block <- met & (ut$U_m > cur_m) & (sweep(ut$U_w, 2, cur_w, ">"))
  block[cbind(couples$man, couples$woman)] <- FALSE
  sum(block)
}
