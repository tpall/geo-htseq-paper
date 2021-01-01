
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
  
}

#+ libs
library(stats) # masks filter
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(gt)
library(extrafont)
library(cowplot)
library(patchwork)
library(viridis)
library(brms)
library(rstan)
library(tidybayes)
library(magick)
library(here)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

#+ params
chains <- 4
cores <- chains
refresh = 200
rstan_options(auto_write = TRUE, javascript = FALSE)
if (!dir.exists("results/models")) {
    dir.create("results/models", recursive = TRUE)
}

#+ data
conformity_acc <- read_csv(here("results/conformity_acc.csv"))
pvalues <- read_csv(here("results/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))

#'
#+ fig1
f <- conforms ~ year
family <- bernoulli()
data <- conformity_acc
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           refresh = refresh,
           chains = chains,
           cores = cores,
           file = here("results/models/conforms_year.rds"))
p <- plot(conditional_effects(mod), plot = FALSE)$year
p <- p + 
  geom_smooth(color = "black") +
  labs(x = "Year", y = "Proportion of submissions conforming\nwith GEO submission guidelines") +
  scale_x_continuous(breaks = seq(2006, 2019, by = 2)) +
  scale_y_continuous(limits = c(0, 1))
ggsave(here("figures/figure_1.tiff"), plot = p, height = 6, width = 10, dpi = 300, units = "cm")

#' 
#+ fig2
parsed_suppfiles <- read_csv(here("results/data/parsed_suppfiles.csv")) %>% 
  filter(str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())


pvalues_acc <- read_csv(here("results/pvalues_acc.csv"))
suppfiles_sample <- parsed_suppfiles %>% 
  inner_join(pvalues_acc)


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
  write_lines(here("results/files_to_check_pi0.txt"))


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

hist_examples <- list(
  Class = c("conservative", "anti-conservative", "uniform", "bimodal", "other"), 
  id = c("GSE63555_Gene_expression_cuffdiff_fpkm.txt.gz", "GSE89511_diff.geneFpkm.exon.glm.LogFc.0.exc.Gapdh.RPMI8226_MCL1.txt.gz", "GSE102826_model_fc.txt.gz", "GSE115649_P20WT_P20M_results.csv.gz", "GSE98869_diffexp.tsv.gz"),
  Set = c("p_value", "vec.v.cd86.pvalue", "p.value (mutunc5 - control)", "pvalue", "pvalue")
) %>% 
  as_tibble()


hist_data <- parsed_suppfiles %>% 
  inner_join(hist_examples)


hist_data_plots <- hist_data %>% 
  mutate(hist = str_remove_all(hist, "[:punct:]"),
         hist = str_split(hist, " "),
         hist = map(hist, as.numeric)) %>% 
  select(Accession, id, Class, hist) %>%
  mutate(QC_thr = map_dbl(hist, qc_threshold, fdr = 0.05),
         QC_plot = map2(hist, QC_thr, plot_qc_hist))


class_counts <- suppfiles_sample %>% 
  count(Class, name = "N")

fit <- brm(Class ~ 1, 
           data = suppfiles_sample, 
           family = categorical(), 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("results/models/Class_1.rds"))

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

# tibble_output %>%
#   gtsave(here("figures/figure_2.pdf"))
# knitr::plot_crop(here("figures/figure_2.pdf"))
tibble_output %>% 
  gtsave(here("figures/figure_2.png"), 
         expand = 10)
fig_2 <- image_read(here("figures/figure_2.png"))
image_write(fig_2, path = here("figures/figure_2.tiff"), format = "tiff")

#'
#+ fig3
f <- Class ~ year + (year | de_tool)
family <- categorical()
data <- pvalues_sample
get_prior(f, data, family)
priors <- c(
  set_prior("lkj(3)", class = "cor"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("normal(0, 0.5)", class = "sd", dpar = "muuniform"),
  set_prior("normal(0, 0.5)", class = "sd", dpar = "mubimodal"),
  set_prior("normal(0, 0.5)", class = "sd", dpar = "muconservative"),
  set_prior("normal(0, 0.5)", class = "sd", dpar = "muother")
  )
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = 3, 
           cores = 3, 
           refresh = refresh,
           prior = priors,
           control = list(adapt_delta = 0.99),
           iter = 4000,
           file = here("results/models/Class_year__year_detool_year.rds"))
conditions <- make_conditions(data, vars = "de_tool")
rownames(conditions) <- conditions$de_tool
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions,
                              categorical = TRUE,
                              re_formula = NULL),
          plot = FALSE)
p3a <- p$`year:cats__` + 
  scale_x_continuous(breaks = seq(2009, 2019, by = 2)) +
  facet_wrap(~de_tool, nrow = 1) +
  labs(y = "Proportion",
       x = "Year") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
p3a$layers[[1]]$aes_params$alpha <- 0.2

#'
#'
#+
f <- Class ~ de_tool
family <- categorical()
data <- pvalues_sample %>% filter(year >= 2018) # keep only values from 2018-2019!!!
get_prior(f, data, family)
priors <- c(
  set_prior("normal(0, 1)", class = "b", dpar = "muuniform")
)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = priors,
           chains = chains, 
           cores = chains, 
           refresh = refresh,
           control = list(adapt_delta = 0.99),
           iter = 2000,
           file = here("results/models/Class_detool_2018-19.rds"))
p <- plot(conditional_effects(mod, 
                              categorical = TRUE, 
                              effects = "de_tool",
                              re_formula = NULL), 
          plot = FALSE)
p3b <- p$`de_tool:cats__` + 
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.position = "bottom")
p3b$layers[[1]]$aes_params$size <- 1
p3 <- p3a / p3b + plot_annotation(tag_levels = "A") +  plot_layout(guides = 'auto')
ggsave(here("figures/figure_3.pdf"), plot = p3, width = 18, height = 12, units = "cm", dpi = 300)

#'
#'
#+
f <- pi0 ~ de_tool
family <- student()
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("results/models/pi0_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool",
                              re_formula = NULL),
          plot = FALSE)
p4b <- p$de_tool +
  labs(y = expression(pi * 0)) +
  theme(axis.title.x = element_blank())
p4b$layers[[1]]$aes_params$size <- 1

#+
f <- pi0 ~ year + (year | de_tool)
prior <- c(prior(normal(0, 0.5), class = b), prior(lkj(3), class = cor), prior(student_t(5, 0, 0.5), class = sigma))
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           prior = prior,
           iter = 3000,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
           file = here("results/models/pi0_year__year_detool.rds"))

conditions <- make_conditions(data, vars = "de_tool")
rownames(conditions) <- conditions$de_tool
p <- plot(conditional_effects(mod, 
                              effects = "year",
                              conditions = conditions,
                              re_formula = NULL),
          plot = FALSE)
p4c <- p$year + 
  facet_wrap(~ de_tool) +
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2009, 2019, 2)) +
  labs(y = expression(pi * 0))
p4c$layers[[1]]$aes_params$alpha <- 0.2

p4a <- pvalues_sample %>% 
  filter(Class == "anti-conservative") %>% 
  ggplot() + 
  geom_histogram(aes(pi0), color = "white", binwidth = 0.1) +
  labs(x = expression(pi * 0), y = "Count")

#+ Fig4, fig.cap=""
p4 <- (p4a + p4b) / p4c + 
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "A")
ggsave(here("figures/figure_4.pdf"), plot = p4, width = 18, height = 12, units = "cm", dpi = 300)
