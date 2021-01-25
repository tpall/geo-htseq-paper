
library("readr")
library("dplyr")
library("ggplot2")
library("here")
library(cowplot)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

data <- read_csv(here("results/data/de_simulation_results.csv"))

qc_threshold <- function(n_pvals, bins, fdr) {
  qbinom(1 - 1 / bins * fdr, n_pvals, 1 / bins)
}

qc_thr <- qc_threshold(n_pvals = 20000, bins = 40, fdr = 0.05)

p <- data %>% 
  mutate(
    reps = factor(reps, levels = c("N=3", "N=6", "N=10")),
    n_eff = factor(n_eff, levels = c("Effects=100", "Effects=200", "Effects=400", "Effects=800"))
    ) %>% 
  ggplot() +
  geom_histogram(aes(pvalue), binwidth = 1/40, center = 1/80) +
  geom_hline(yintercept = qc_thr, color = "red") +
  facet_grid(reps~n_eff) +
  scale_x_continuous(breaks = c(0, 0.5, 1))
ggsave(here("figures/figure_2_figure_supplement_1.tiff"), plot = p, width = 12, height = 8, units = "cm", dpi = 300)

