# Tests for the core library. Run from the repository root:
#   Rscript tests/test_core.R
for (f in c("population", "matching", "estimate", "run_helpers")) source(file.path("R", paste0(f, ".R")))
set.seed(1)
fails <- 0
check <- function(ok, msg) { cat(if (ok) "PASS " else "FAIL ", msg, "\n"); if (!ok) fails <<- fails + 1 }
b <- c(same_group = 1, edu_diff = -1)

# Gale-Shapley: one-to-one, individually rational, stable, deterministic.
for (u0 in c(-Inf, 0, 1.5)) for (r in 1:3) {
  pop <- make_population(80); ut <- make_utilities(pop, b, b, u0); cp <- gale_shapley(ut)
  tag <- sprintf("u0 = %s, market %d:", u0, r)
  check(!anyDuplicated(cp$man) && !anyDuplicated(cp$woman), paste(tag, "one-to-one"))
  check(all(ut$U_m[cbind(cp$man, cp$woman)] > ut$out_m[cp$man]) &&
        all(ut$U_w[cbind(cp$man, cp$woman)] > ut$out_w[cp$woman]),
        paste(tag, "both prefer the union to staying single"))
  check(count_blocking_pairs(ut, cp) == 0, paste(tag, "no blocking pairs"))
  check(identical(cp, gale_shapley(ut)), paste(tag, "same utilities give the same matching"))
}
pop <- make_population(80); ut <- make_utilities(pop, b, b, -Inf)
check(nrow(gale_shapley(ut)) == 80, "with no outside option everyone pairs")

# One-sided choice picks each person's best option.
pop <- make_population(50); ut <- make_utilities(pop, b, b, -Inf)
ch <- one_sided(ut$U_m, ut$out_m)
check(all(ch == apply(ut$U_m, 1, which.max)), "one-sided choice is each chooser's best option")

# Calibration hits the target pairing rate.
pop <- make_population(300)
mk <- function(u) make_utilities(pop, b, b, u)
u0 <- calibrate_u0(mk, 0.85)
set.seed(9001); rate <- nrow(gale_shapley(mk(u0))) / 300
check(abs(rate - 0.85) < 0.03, sprintf("calibrated u0 = %.2f gives pairing %.3f (target 0.85)", u0, rate))

# Choice sets: partner first, J alternatives, no duplicates within a set.
pop <- make_population(200); cs <- build_choice_sets(1:20, 1:20, pop$men, pop$women, 30)
check(nrow(cs) == 600 && all(tapply(cs$chosen, cs$chooser, sum) == 1), "one chosen row per choice set")

# The predicted share equals the observed share at the estimates.
pop <- make_population(1000); ut <- make_utilities(pop, b, b, -Inf)
ch <- one_sided(ut$U_m, ut$out_m); cs <- build_choice_sets(1:1000, ch, pop$men, pop$women, 100)
f <- fit_clogit(cs, c("same_group", "alt_minority", "edu_diff"))
check(abs(predicted_share(cs, setNames(f$est, f$term)) - mean(cs$same_group[cs$chosen == 1])) < 1e-6,
      "predicted share equals observed share at the estimates")

cat(if (fails == 0) "\nAll tests passed.\n" else sprintf("\n%d FAILED\n", fails))
quit(status = as.integer(fails > 0))
