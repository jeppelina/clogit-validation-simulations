# =============================================================================
# run_helpers.R: the two settings and one simulated market
# =============================================================================
#   one_sided   each person picks freely; no outside option
#   two_sided   Gale-Shapley; the value of staying single is set so that
#               85% of people pair up
# =============================================================================

SETTINGS <- c("one_sided", "two_sided")
TARGET_PAIRING <- 0.85

#' u0 for a setting. make_ut(u0) builds a market from fixed random numbers.
setting_u0 <- function(setting, make_ut)
  if (setting == "one_sided") -Inf else calibrate_u0(make_ut, TARGET_PAIRING)

#' Form unions. Returns, for each sex, the egos who are partnered and their
#' partners, plus the couples and the pairing rate. Under one-sided choice
#' each sex's choices are its own (a person can be chosen by many).
form_unions <- function(pop, ut, setting) {
  N <- nrow(pop$men)
  if (setting == "one_sided") {
    ch_m <- one_sided(ut$U_m, ut$out_m); ch_w <- one_sided(t(ut$U_w), ut$out_w)
    men <- which(!is.na(ch_m)); wom <- which(!is.na(ch_w))
    return(list(
      sides = list(men = list(ego = men, alt = ch_m[men], edf = pop$men, adf = pop$women),
                   women = list(ego = wom, alt = ch_w[wom], edf = pop$women, adf = pop$men)),
      couples = data.frame(man = men, woman = ch_m[men]), pairing = length(men) / N))
  }
  cp <- gale_shapley(ut)
  list(sides = list(men = list(ego = cp$man, alt = cp$woman, edf = pop$men, adf = pop$women),
                    women = list(ego = cp$woman, alt = cp$man, edf = pop$women, adf = pop$men)),
       couples = cp, pairing = nrow(cp) / N)
}

#' Start a PSOCK cluster with the core library loaded. (Forked workers can
#' crash survival::clogit on macOS.)
start_cluster <- function(n = max(1, parallel::detectCores() - 2)) {
  cl <- parallel::makeCluster(n)
  invisible(parallel::clusterEvalQ(cl, { for (f in c("population", "matching", "estimate", "run_helpers"))
    source(file.path("R", paste0(f, ".R"))); NULL }))
  cl
}
