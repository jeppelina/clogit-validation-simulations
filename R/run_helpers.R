# =============================================================================
# run_helpers.R -- one simulated market under a named setting
# =============================================================================
# SETTINGS used across Modules B to D (validated in Module A):
#   oracle  one-sided choice, no outside option (recovers own tastes)
#   gs      Gale-Shapley, men proposing, everyone meets everyone,
#           outside option calibrated to 85% pairing
#   gs_w    as gs, women proposing
#   seq     sequential search: budget 50 over 10 rounds, 85% pairing
# =============================================================================

SETTINGS <- list(
  oracle = list(mech = "oracle"),
  gs     = list(mech = "da", budget = Inf, rounds = 1, proposer = "men"),
  gs_w   = list(mech = "da", budget = Inf, rounds = 1, proposer = "women"),
  seq    = list(mech = "da", budget = 50, rounds = 10, proposer = "men")
)
TARGET_PAIR <- 0.85

# Calibrate u0 for a setting, given a function u -> utilities object that
# builds a fresh market (population + utilities) from fixed random numbers.
setting_u0 <- function(st, make_ut) {
  if (st$mech == "oracle") return(-Inf)
  calibrate_u0(make_ut, st$budget, st$rounds, TARGET_PAIR)
}

# Run partnering. Returns list(sides = list(men, women), pair_rate, mean_met).
# Each side: ego ids, partner ids, ego and alter data frames.
form_unions <- function(pop, ut, st) {
  N <- nrow(pop$men)
  if (st$mech == "oracle") {
    ch_m <- oracle_choice(ut$U_m, ut$out_m); ch_w <- oracle_choice(t(ut$U_w), ut$out_w)
    men <- which(!is.na(ch_m)); wom <- which(!is.na(ch_w))
    return(list(
      sides = list(men = list(ego = men, alt = ch_m[men], edf = pop$men, adf = pop$women),
                   women = list(ego = wom, alt = ch_w[wom], edf = pop$women, adf = pop$men)),
      couples = data.frame(man = men, woman = ch_m[men]),
      pair_rate = length(men) / N, mean_met = N))
  }
  sm <- search_market(ut, st$budget, st$rounds, st$proposer)
  cp <- sm$couples
  list(sides = list(men = list(ego = cp$man, alt = cp$woman, edf = pop$men, adf = pop$women),
                    women = list(ego = cp$woman, alt = cp$man, edf = pop$women, adf = pop$men)),
       couples = cp, pair_rate = nrow(cp) / N, mean_met = mean(rowSums(sm$met)))
}

# Fit one clogit per ego sex. Returns long data.frame (ego_sex, term, est, se).
fit_both_sides <- function(un, terms, J = 100) {
  do.call(rbind, lapply(names(un$sides), function(sx) {
    s <- un$sides[[sx]]
    cs <- build_choice_sets(s$ego, s$alt, s$edf, s$adf, J)
    f <- tryCatch(fit_clogit(cs, terms),
                  error = function(e) data.frame(term = terms, est = NA, se = NA))
    f$ego_sex <- sx
    f
  }))
}
