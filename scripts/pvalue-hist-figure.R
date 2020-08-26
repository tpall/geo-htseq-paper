library(tidyverse)
library(gt)


parsed_suppfiles <- read_csv("data/parsed_suppfiles.csv") %>% 
  filter(str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())


files_to_check <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  filter(near(pi0, 1)) %>% 
  select(id, pi0, FDR_pval, Set) %>% 
  pull(id) %>% 
  str_split(" from ") %>% 
  flatten_chr()


files_to_check[str_detect(files_to_check, "^GSE")] %>% 
  unique() %>% 
  str_c("output/suppl/", .) %>% 
  write_lines("output/files_to_check_pi0.txt")


qc_threshold <- function(x, fdr) {
  bins <- length(x)
  qbinom(1 - 1 / bins * fdr, sum(x), 1 / bins)
}


plot_qc_hist <- function(counts, t) {
  ggplot() +
    geom_col(aes(x = seq(0, 1, length.out = length(counts)), y = counts)) +
    geom_hline(yintercept = t, color = "red", size = 3) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}


set.seed(15)
suppfiles_sample <- parsed_suppfiles %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  ungroup()


hist_data <- suppfiles_sample %>%
  filter((Class %in% c("bimodal", "conservative", "other", "uniform")) | (Class %in% c("anti-conservative") & pi0 >= 0.8)) %>% 
  group_by(Class) %>% 
  sample_n(1)


hist_data_plots <- hist_data %>% 
  mutate(hist = str_remove_all(hist, "[:punct:]"),
         hist = str_split(hist, " "),
         hist = map(hist, as.numeric)) %>% 
  select(Accession, id, Class, hist) %>%
  mutate(QC_thr = map_dbl(hist, qc_threshold, fdr = 0.05),
         QC_plot = map2(hist, QC_thr, plot_qc_hist))


class_counts <- suppfiles_sample %>% 
  count(Class, name = "N")


library(brms)
library(tidybayes)
fit <- brm(Class ~ 1, data = suppfiles_sample, family = categorical())
pe <- posterior_epred(fit)
classes_props <- pe[1:4000, 1, 1:5] %>% 
  as_tibble() %>% 
  map(mean_hdi) %>% 
  bind_rows(.id = "Class") %>% 
  mutate_at(vars(y, ymin, ymax), signif, digits = 2) %>% 
  mutate(`Fraction [95% CI]` = str_c(y, " [", ymin, "; ", ymax, "]"),
         Example = NA) %>% 
  select(Class, `Fraction [95% CI]`, Example)


plot_data <- hist_data_plots %>% 
  select(Class, QC_plot) %>% 
  left_join(class_counts) %>% 
  left_join(
    classes_props
    ) %>% 
  arrange(desc(N)) %>% 
  select(Class, N, `Fraction [95% CI]`, Example, QC_plot) %>% 
  ungroup()


tibble_output <- plot_data %>%
  select(-QC_plot) %>% 
  gt() %>%
  text_transform(
    locations = cells_body(vars(Example)),
    fn = function(x) {
      map(plot_data$QC_plot, ggplot_image, height = px(30), aspect_ratio = 1.5)
    }
    )


tibble_output %>% 
  gtsave("output/hist_examples.png", 
         expand = 10)
