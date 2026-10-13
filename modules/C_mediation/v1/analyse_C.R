# =============================================================================
# Module C: summarise KHB shares. Run from the repo root after run_C.R.
# =============================================================================
OUT <- "modules/C_mediation/results"
r <- read.csv(file.path(OUT, "C.csv"))
mcse <- function(x) sd(x) / sqrt(length(x))
s <- do.call(rbind, lapply(split(r, r[, c("cell", "setting", "spec")], drop = TRUE), function(d)
  data.frame(d[1, c("cell", "setting", "spec")], share = mean(d$share), share_mcse = mcse(d$share),
             b_reduced = mean(d$b_reduced), b_full = mean(d$b_full), b_mediator = mean(d$b_mediator),
             cor_med_focal = mean(d$cor_med_focal), pair_rate = mean(d$pair_rate))))
rownames(s) <- NULL
# reference: oracle share in the same cell (researcher spec), and for C3 the
# trait-controlled share under the same setting
ref <- s[s$setting == "oracle" & s$spec == "as_researcher", c("cell", "share")]
names(ref)[2] <- "oracle_share"
s <- merge(s, ref, by = "cell")
s$gap_vs_oracle <- s$share - s$oracle_share
tc <- s[s$spec == "trait_controlled", c("cell", "setting", "share")]
names(tc)[3] <- "trait_controlled_share"
s <- merge(s, tc, by = c("cell", "setting"), all.x = TRUE)
s$inflation <- ifelse(s$spec == "as_researcher", s$share - s$trait_controlled_share, NA)
s <- s[order(s$cell, s$spec, match(s$setting, c("oracle", "gs", "seq"))), ]
write.csv(s, file.path(OUT, "C_summary.csv"), row.names = FALSE)
options(width = 180)
print(s[, c("cell", "setting", "spec", "share", "share_mcse", "gap_vs_oracle", "inflation",
            "b_reduced", "b_full", "b_mediator", "cor_med_focal")], digits = 3, row.names = FALSE)
