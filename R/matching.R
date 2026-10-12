# =============================================================================
# matching.R -- how partners are formed
# =============================================================================
# oracle_choice(): each ego picks the best option in the whole opposite-sex
#   pool, or stays single. No capacity limits, the other side has no say.
#   This is exactly the process conditional logit assumes. Run per sex.
#
# search_market(): sequential search with deferred acceptance (the default
#   two-sided mechanism).
#   - Tastes are fixed: utilities and outside options are drawn once.
#   - Each round, every person who is currently single meets new people of
#     the other sex at random. Meetings are mutual and accumulate.
#   - After each round the market clears by deferred acceptance among
#     everyone who has met so far. Couples are provisional until the last
#     round, so a later round can break them up.
#   - Rounds stop when everyone is matched or the rounds are used up.
#   - The final matching is stable with respect to the final meeting network:
#     no two people who met would both rather be with each other.
#   The search budget is the expected number of people a person meets if
#   single throughout. With one round this is a one-shot meeting network.
#   With budget >= N everyone meets everyone, which is textbook Gale-Shapley.
# =============================================================================

oracle_choice <- function(U_ego, out_ego) {
  # U_ego: [ego, alter] utilities. Returns the chosen alter per ego, NA = single.
  best <- max.col(U_ego, ties.method = "first")
  best_u <- U_ego[cbind(seq_len(nrow(U_ego)), best)]
  ifelse(best_u > out_ego, best, NA_integer_)
}

# Deferred acceptance on a fixed meeting network.
#   U_prop[p, r]: proposer p's utility of receiver r
#   U_recv[p, r]: receiver r's utility of proposer p
#   met[p, r]:    TRUE if p and r have met
# Returns the receiver each proposer ends up with (NA = single).
deferred_acceptance <- function(U_prop, U_recv, out_prop, out_recv, met) {
  n_p <- nrow(U_prop); n_r <- ncol(U_prop)
  # Each proposer's list: met receivers preferred to staying single, best first.
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

# Clear the market on a meeting network, either side proposing.
# Returns a data.frame of couples (man, woman).
clear_market <- function(ut, met, proposer = "men") {
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

# Each person in `who` meets new people they have not met, drawn at random
# from the other side. n_new is a per-person count: either one number or one
# per person in `who`. side = "men" draws women for men; "women" the reverse.
add_meetings <- function(met, who, n_new, side) {
  if (length(who) == 0) return(met)
  n_new <- rep_len(n_new, length(who))
  for (j in seq_along(who)) {
    i <- who[j]; k <- n_new[j]
    if (k <= 0) next
    if (side == "men") {
      cand <- which(!met[i, ])
      if (length(cand)) met[i, cand[sample.int(length(cand), min(k, length(cand)))]] <- TRUE
    } else {
      cand <- which(!met[, i])
      if (length(cand)) met[cand[sample.int(length(cand), min(k, length(cand)))], i] <- TRUE
    }
  }
  met
}

# Sequential search. budget: expected meetings per person if single
# throughout. rounds: number of search rounds (1 = one-shot). Each round,
# each single man and each single woman draws budget / (2 * rounds) new
# contacts, so a person single throughout meets about `budget` people in all
# (their own draws plus draws that land on them).
search_market <- function(ut, budget, rounds = 10, proposer = "men") {
  N_m <- nrow(ut$U_m); N_w <- ncol(ut$U_m)
  if (budget >= min(N_m, N_w)) {
    met <- matrix(TRUE, N_m, N_w)
    return(list(couples = clear_market(ut, met, proposer), met = met,
                 rounds_used = 1L))
  }
  per_draw <- budget / (2 * rounds)
  met <- matrix(FALSE, N_m, N_w)
  couples <- data.frame(man = integer(0), woman = integer(0))
  rounds_used <- 0L
  for (r in seq_len(rounds)) {
    single_m <- setdiff(seq_len(N_m), couples$man)
    single_w <- setdiff(seq_len(N_w), couples$woman)
    if (length(single_m) == 0 && length(single_w) == 0) break
    # Fractional draws are resolved per person, so the expected count is
    # exact and no round is empty for everyone at once.
    draws <- function(n) floor(per_draw) + (runif(n) < per_draw - floor(per_draw))
    met <- add_meetings(met, single_m, draws(length(single_m)), "men")
    met <- add_meetings(met, single_w, draws(length(single_w)), "women")
    couples <- clear_market(ut, met, proposer)
    rounds_used <- r
  }
  list(couples = couples, met = met, rounds_used = rounds_used)
}

# Find the outside-option level u0 that gives a target pairing rate.
# make_ut(u0) must return a utilities object built with that u0 and FIXED
# random numbers (same seed each call), so the rate is monotone in u0.
calibrate_u0 <- function(make_ut, budget, rounds, target = 0.85,
                         lower = -5, upper = 15, iters = 14, n_markets = 2) {
  rate_at <- function(u0) {
    mean(vapply(seq_len(n_markets), function(s) {
      set.seed(9000 + s)
      ut <- make_ut(u0)
      sm <- search_market(ut, budget, rounds)
      nrow(sm$couples) / nrow(ut$U_m)
    }, numeric(1)))
  }
  for (i in seq_len(iters)) {
    mid <- (lower + upper) / 2
    if (rate_at(mid) > target) lower <- mid else upper <- mid
  }
  (lower + upper) / 2
}

# Stability check: number of blocking pairs among people who met
# (both strictly prefer each other to their current outcome).
count_blocking_pairs <- function(ut, couples, met) {
  cur_m <- ut$out_m; cur_w <- ut$out_w
  cur_m[couples$man] <- ut$U_m[cbind(couples$man, couples$woman)]
  cur_w[couples$woman] <- ut$U_w[cbind(couples$man, couples$woman)]
  block <- met & (ut$U_m > cur_m) & (sweep(ut$U_w, 2, cur_w, ">"))
  block[cbind(couples$man, couples$woman)] <- FALSE
  sum(block)
}

# --- Legacy: one-shot network where only men draw k women. Used by the
# --- 2026-10-12 A0 pilot (pilot_A0.R); kept so that pilot reproduces.
meeting_network <- function(N_m, N_w, k) {
  if (k >= N_w) return(matrix(TRUE, N_m, N_w))
  M <- matrix(FALSE, N_m, N_w)
  for (m in seq_len(N_m)) M[m, sample.int(N_w, k)] <- TRUE
  M
}
match_market <- function(ut, k, proposer = c("men", "women")) {
  proposer <- match.arg(proposer)
  clear_market(ut, meeting_network(nrow(ut$U_m), ncol(ut$U_m), k), proposer)
}
