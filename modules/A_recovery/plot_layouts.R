# =============================================================================
# Module A (v3): plot the residential layouts used in test 2 and report how
# segregated each is. Run from the repo root.
# For each layout: one market (1,000 men, 1,000 women, 80/20), and two
# measures over the opposite sex, i.e. the partner market:
#   minority exposure: mean share of the minority among each minority
#     person's 50 nearest people of the other sex (the 20% line = no
#     segregation)
#   dissimilarity index on a 10 x 10 grid (0 = even, 1 = complete)
# =============================================================================
library(ggplot2)
source("R/population.R")
set.seed(20261014)
layouts <- list(c("mixed", 0), c("sides", 0.5), c("sides", 0.9), c("enclaves", 0.5), c("enclaves", 0.9))
labels <- c("Mixed", "Two sides, moderate", "Two sides, strong", "Enclaves, moderate", "Enclaves, strong")
pts <- list(); stats <- list()
for (i in seq_along(layouts)) {
  pop <- make_population(1000, c(0.8, 0.2), layout = layouts[[i]][1], s = as.numeric(layouts[[i]][2]))
  m <- pop$men; w <- pop$women
  d <- sqrt(outer(m$x, w$x, "-")^2 + outer(m$y, w$y, "-")^2)
  mi <- which(m$group == 2)
  expo <- mean(apply(d[mi, , drop = FALSE], 1, function(r) mean(w$group[order(r)[1:50]] == 2)))
  cell <- function(p) paste(pmin(floor(p$x * 10), 9), pmin(floor(p$y * 10), 9))
  all <- rbind(m, w); cl <- factor(cell(all))
  a <- table(cl[all$group == 1]) / sum(all$group == 1); b <- table(cl[all$group == 2]) / sum(all$group == 2)
  dis <- 0.5 * sum(abs(a - b))
  stats[[i]] <- data.frame(layout = labels[i], minority_exposure = expo, dissimilarity = dis)
  pts[[i]] <- data.frame(layout = labels[i], x = m$x, y = m$y,
                         group = ifelse(m$group == 1, "majority (80%)", "minority (20%)"))
}
stats <- do.call(rbind, stats); print(stats, digits = 2, row.names = FALSE)
pts <- do.call(rbind, pts)
lab <- setNames(sprintf("%s\nexposure %.2f, dissimilarity %.2f", stats$layout, stats$minority_exposure,
                        stats$dissimilarity), stats$layout)
pts$layout <- factor(pts$layout, labels, lab[labels])
p <- ggplot(pts, aes(x, y, colour = group)) +
  geom_point(size = 0.6, alpha = 0.7) + coord_equal() +
  facet_wrap(~ layout, nrow = 1) +
  scale_colour_manual(values = c("grey60", "#D55E00")) +
  labs(x = NULL, y = NULL, colour = NULL, title = "Residential layouts in test 2 (men shown; women are placed the same way)",
       subtitle = "Exposure: minority share among a minority person's 50 nearest people of the other sex (0.20 = none)") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", axis.text = element_blank(),
                                       panel.grid = element_blank(), panel.border = element_rect(fill = NA, colour = "grey80"))
dir.create("modules/A_recovery/results", showWarnings = FALSE)
ggsave("modules/A_recovery/results/layouts.pdf", p, width = 13, height = 3.6)
write.csv(stats, "modules/A_recovery/results/layouts_stats.csv", row.names = FALSE)
