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

cat(if (fails == 0) "\nAll tests passed.\n" else sprintf("\n%d FAILED\n", fails))
quit(status = as.integer(fails > 0))
