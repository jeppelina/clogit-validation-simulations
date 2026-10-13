# =============================================================================
# Study A: draw the residential layouts used in test 2 and measure how
# segregated each is. Run from the repository root:
#   Rscript A_recovery/plot_layouts.R
# Measures, over the other sex (the partner market):
#   exposure:      minority share among a minority person's 50 nearest
#                  people (0.20 = no segregation)
#   dissimilarity: on a 10 x 10 grid (0 = even, 1 = complete separation)
# =============================================================================
library(ggplot2)
# Save a figure as PDF and as PNG (the PNG is shown in the README).
save_fig <- function(file, plot, width, height) {
  ggsave(file, plot, width = width, height = height)
  ggsave(sub("\\.pdf$", ".png", file), plot, width = width, height = height, dpi = 150, bg = "white")
}
source("R/population.R")
set.seed(1)
layouts <- list(c("mixed", 0), c("sides", 0.5), c("sides", 0.9), c("enclaves", 0.5), c("enclaves", 0.9))
labels <- c("Mixed", "Two sides, moderate", "Two sides, strong", "Enclaves, moderate", "Enclaves, strong")
pts <- list(); stats <- list()
for (i in seq_along(layouts)) {
  pop <- make_population(1000, c(0.8, 0.2), layout = layouts[[i]][1], s = as.numeric(layouts[[i]][2]))
  m <- pop$men; w <- pop$women
  d <- sqrt(outer(m$x, w$x, "-")^2 + outer(m$y, w$y, "-")^2)
  expo <- mean(apply(d[m$group == 2, , drop = FALSE], 1, function(r) mean(w$group[order(r)[1:50]] == 2)))
  all <- rbind(m, w)
  cl <- factor(paste(pmin(floor(all$x * 10), 9), pmin(floor(all$y * 10), 9)))
  a <- table(cl[all$group == 1]) / sum(all$group == 1); b <- table(cl[all$group == 2]) / sum(all$group == 2)
  stats[[i]] <- data.frame(layout = labels[i], exposure = expo, dissimilarity = 0.5 * sum(abs(a - b)))
  pts[[i]] <- data.frame(layout = labels[i], x = m$x, y = m$y,
                         group = ifelse(m$group == 1, "majority (80%)", "minority (20%)"))
}
stats <- do.call(rbind, stats); print(stats, digits = 2, row.names = FALSE)
pts <- do.call(rbind, pts)
lab <- setNames(sprintf("%s\nexposure %.2f, dissimilarity %.2f", stats$layout, stats$exposure, stats$dissimilarity), stats$layout)
pts$layout <- factor(pts$layout, labels, lab[labels])
p <- ggplot(pts, aes(x, y, colour = group)) + geom_point(size = 0.6, alpha = 0.7) + coord_equal() +
  facet_wrap(~ layout, nrow = 1) + scale_colour_manual(values = c("grey60", "#D55E00")) +
  labs(x = NULL, y = NULL, colour = NULL, title = "Residential layouts (men shown; women are placed the same way)",
       subtitle = "Dissimilarity above 0 in the mixed layout is sampling noise (about 20 people per grid cell).") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", axis.text = element_blank(),
    panel.grid = element_blank(), panel.border = element_rect(fill = NA, colour = "grey80"))
dir.create("A_recovery/results", showWarnings = FALSE)
save_fig("A_recovery/results/layouts.pdf", p, width = 13, height = 3.6)
write.csv(stats, "A_recovery/results/layouts.csv", row.names = FALSE)
