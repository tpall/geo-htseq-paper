library(tidyverse)

parsed_suppfiles <- read_csv("data/parsed_suppfiles.csv") %>% 
  filter(str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())

qc_threshold <- function(x, fdr) {
  bins <- length(x)
  qbinom(1 - 1 / bins * fdr, sum(x), 1 / bins)
}

plot_qc_hist <- function(counts, t) {
  ggplot() +
    geom_col(aes(x = seq(0, 1, length.out = length(counts)), y = counts)) +
    geom_hline(yintercept = t, color = "red") +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

set.seed(12)
hist_data <- parsed_suppfiles %>% 
  filter((Class %in% c("bimodal", "conservative", "other", "uniform")) | (Class %in% c("anti-conservative") & pi0 >= 0.8)) %>% 
  group_by(Class) %>% 
  sample_n(1)

hist_data_plots <- hist_data %>% 
  mutate(hist = str_split(hist, "[^0-9]+"),
         hist = map(hist, as.numeric),
         hist = map(hist, na.omit)) %>% 
  select(Accession, id, Class, hist) %>%
  mutate(QC_thr = map_dbl(hist, qc_threshold, fdr = 0.05),
         QC_plot = map2(hist, QC_thr, plot_qc_hist))

hist_data_plots$QC_plot

tibble_output <- tibble(
  text = hist_data_plots$Class,
  ggplot = NA,
  .rows = length(hist_data_plots$Class)) %>%
  gt() %>%
  text_transform(
    locations = cells_body(vars(ggplot)),
    fn = function(x) {
      map(hist_data_plots$QC_plot, ggplot_image, height = px(100), aspect_ratio = 1.5)
    }
    )

tibble_output
