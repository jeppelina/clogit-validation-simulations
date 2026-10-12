# Tests for R/matching.R. Run from the repo root: Rscript tests/test_matching.R
source("R/population.R"); source("R/matching.R")
set.seed(1)
fails <- 0
check <- function(ok, msg) { cat(if (ok) "PASS " else "FAIL ", msg, "\n"); if (!ok) fails <<- fails + 1 }

beta <- c(same_group = 1, edu_diff = -1)
for (k in c(5, 20, 60)) for (u0 in c(-Inf, 1)) for (prop in c("men", "women")) {
  pop <- make_population(60)
  ut <- make_utilities(pop, beta, beta, u0)
  met <- meeting_network(60, 60, k)
  # match_market draws its own network; rebuild the call so we can check
  # stability against the same network
  if (prop == "men") {
    pm <- deferred_acceptance(ut$U_m, ut$U_w, ut$out_m, ut$out_w, met)
    cp <- data.frame(man = which(!is.na(pm)), woman = pm[!is.na(pm)])
  } else {
    pw <- deferred_acceptance(t(ut$U_w), t(ut$U_m), ut$out_w, ut$out_m, t(met))
    cp <- data.frame(man = pw[!is.na(pw)], woman = which(!is.na(pw)))
  }
  tag <- sprintf("k=%d u0=%s %s-proposing", k, u0, prop)
  check(!anyDuplicated(cp$man) && !anyDuplicated(cp$woman), paste(tag, "one-to-one"))
  check(all(met[cbind(cp$man, cp$woman)]), paste(tag, "couples only among people who met"))
  check(all(ut$U_m[cbind(cp$man, cp$woman)] > ut$out_m[cp$man]) &&
        all(ut$U_w[cbind(cp$man, cp$woman)] > ut$out_w[cp$woman]),
        paste(tag, "both prefer the union to staying single"))
  check(count_blocking_pairs(ut, cp, met) == 0, paste(tag, "no blocking pairs"))
}

# Fixed tastes: same utilities in, same matching out.
pop <- make_population(80); ut <- make_utilities(pop, beta, beta, 0)
met <- meeting_network(80, 80, 80)
a <- deferred_acceptance(ut$U_m, ut$U_w, ut$out_m, ut$out_w, met)
b <- deferred_acceptance(ut$U_m, ut$U_w, ut$out_m, ut$out_w, met)
check(identical(a, b), "deterministic given utilities (tastes are fixed)")

# Sequential search: final matching is stable on the final network, couples
# only among people who met, individually rational; budget >= N equals
# textbook Gale-Shapley on the complete network.
for (budget in c(4, 12, 30)) for (rounds in c(1, 5)) for (u0 in c(-Inf, 1)) {
  pop <- make_population(60); ut <- make_utilities(pop, beta, beta, u0)
  sm <- search_market(ut, budget, rounds); cp <- sm$couples
  tag <- sprintf("search budget=%d rounds=%d u0=%s", budget, rounds, u0)
  check(!anyDuplicated(cp$man) && !anyDuplicated(cp$woman), paste(tag, "one-to-one"))
  check(all(sm$met[cbind(cp$man, cp$woman)]), paste(tag, "couples only among people who met"))
  check(count_blocking_pairs(ut, cp, sm$met) == 0, paste(tag, "no blocking pairs on final network"))
  check(sm$rounds_used >= 1 && sm$rounds_used <= rounds, paste(tag, "rounds within budget"))
}
pop <- make_population(60); ut <- make_utilities(pop, beta, beta, 0)
gs <- clear_market(ut, matrix(TRUE, 60, 60))
sm <- search_market(ut, 60, 10)
check(identical(gs, sm$couples), "budget >= N reproduces textbook Gale-Shapley")

# Expected meetings: a person single throughout meets about `budget` people.
pop <- make_population(400); ut <- make_utilities(pop, beta, beta, 100)  # nobody pairs
sm <- search_market(ut, 40, 10)
check(abs(mean(rowSums(sm$met)) - 40) < 2, sprintf("mean meetings when all single = %.1f (target 40)", mean(rowSums(sm$met))))

# Calibration hits its target pairing rate.
pop <- make_population(300)
mk <- function(u0) make_utilities(pop, beta, beta, u0)
u0 <- calibrate_u0(mk, budget = 15, rounds = 10, target = 0.85)
set.seed(9001); r1 <- nrow(search_market(mk(u0), 15, 10)$couples) / 300
check(abs(r1 - 0.85) < 0.03, sprintf("calibrated u0 = %.2f gives pairing rate %.3f (target 0.85)", u0, r1))

cat(if (fails == 0) "\nAll tests passed.\n" else sprintf("\n%d FAILED\n", fails))
quit(status = as.integer(fails > 0))
