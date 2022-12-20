#' ---
#' title: Supporting Information
#' author: ""
#' date: ""
#' output:
#'    bookdown::pdf_document2:
#'       toc: FALSE
#' urlcolor: blue
#' header-includes:
#' - \usepackage[font=footnotesize,labelfont=bf]{caption}
#' - \DeclareCaptionLabelFormat{supplement}{S#2 Fig.}
#' - \captionsetup{labelformat=supplement,labelsep=quad}
#' bibliography: main/references.bib
#' csl: main/plos-biology.csl
#' ---


#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, message = FALSE, comment = FALSE, warning = FALSE)
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}


#+ libs
library(stats) # masks filter
library(readr)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(stringr)
library(cowplot)
library(patchwork)
library(brms)
library(rstan)
library(tidybayes)
library(forcats)
library(here)
library(glue)
library(modelr)
old <- theme_set(theme_cowplot(font_size = 10, font_family = "Helvetica"))

#'
#+ params
is_ci <- function() {
  "CI" %in% Sys.getenv()
}
chains <- ifelse(is_ci(), 1, 4)
cores <- chains
refresh <- 0
rstan_options(auto_write = TRUE, javascript = FALSE)
if (!dir.exists(here("results/models"))) {
  message("Creating results/models dir..")
  dir.create(here("results/models"), recursive = TRUE)
}

#+ data
parse_detool <- . %>% 
  mutate(
    de_tool = case_when(
      is.na(de_tool) ~ "unknown",
      de_tool %in% c("biojupie", "exdega", "nascent rna seq") ~ "unknown",
      de_tool == "cufflinks-edger" ~ "edger",
      de_tool == "cufflinks-deseq" ~ "deseq",
      de_tool == "cufflinks-deseq2" ~ "deseq2",
      de_tool %in% c("cufflinks-limma", "edger-limma") ~ "limma",
      de_tool == "clc genomics-edger" ~ "clc genomics",
      TRUE ~ de_tool
    )
  )

pvalues <- read_csv(here("results/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform) %>% 
  parse_detool()
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform) %>% 
  parse_detool()
pvalues_filtered <- read_csv(here("results/pvalues_filtered.csv")) %>% 
  rename(de_tool = analysis_platform) %>% 
  parse_detool()
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))
de_simulation_results <- read_csv(here("results/data/de_simulation_results.csv"))

#
#+
spots <- read_csv(here("output/spots.csv")) #Parsing proxy for number of samples
n_libs <- spots %>% 
  rename(Accession = geo_accession) %>% 
  count(Accession, name = "n_libs")
n_pvsets <- pvalues %>% 
  count(Accession, name = "n_pvsets")
n_samples <- n_libs %>% 
  inner_join(n_pvsets) %>% 
  mutate(N = n_libs / (n_pvsets + 1))

#'
#'
#+ pi0detoolsample, include=FALSE
priors <- set_prior("normal(0, 2)", class="b")
data <- pvalues_sample %>% 
  drop_na(pi0, de_tool)
mod <- brm(
  formula = pi0 ~ de_tool, 
  data = data, 
  family = Beta(),
  prior = priors,
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/pi0_detool_sample.rds"),
  file_refit = "on_change"
)

#' 
#+ s1fig
get_n_tests <- function(x) {
  x %>% 
    str_remove_all("[\\[\\]]") %>% 
    str_split(",", simplify = TRUE) %>% 
    str_trim() %>% 
    as.numeric() %>% 
    sum(na.rm = TRUE)
}

number_of_pvalues <- pvalues_sample %>%
  select(Accession, id, hist, Class) %>%
  rowwise() %>%
  mutate(
    n_pvalues = map_dbl(hist, get_n_tests),
    anticons = if_else(Class %in% c("anti-conservative", "uniform"), "anti-conservative", "all other classes")
  )

pa <- number_of_pvalues %>% 
  mutate(
    anticons = ifelse(anticons == "anti-conservative", "anti-conservative\nand uniform", anticons)) %>% 
  ggplot() +
  geom_histogram(aes(n_pvalues), binwidth = 0.1) +
  geom_vline(xintercept = median(number_of_pvalues$n_pvalues, na.rm = TRUE), linetype = "dashed") +
  facet_wrap(~anticons, scales = "free_y") +
  scale_x_log10() +
  labs(x = expression(Number~of~p~values~(log[10])), y = "Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

data <- number_of_pvalues %>% 
  mutate(log10_n_pvalues = log10(n_pvalues)) %>% 
  drop_na(log10_n_pvalues, anticons)

mod <- brm(
  log10_n_pvalues ~ anticons,
  data = data,
  family = student(),
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/log10_n_pvalues ~ anticons.rds"),
  file_refit = "on_change"
)

draws <- data %>% 
  data_grid(anticons) %>% 
  add_epred_draws(mod) %>% 
  mutate(anticons = ifelse(anticons == "anti-conservative", "anti-conservative\nand uniform", "all\nother\nclasses"))
pb <- draws %>% 
  ggplot(aes(anticons, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(Number~of~p~values~(log[10])), y = "Count") +
  theme(axis.title.x = element_blank())

p <- pa + pb + plot_layout(widths = c(2/3, 1/3)) + plot_annotation(tag_levels = "A")

#+ include=FALSE
# tiff(here("figures/S1Fig.tiff"), width = 12, height = 8, units = "cm", res = 300)
# p
# dev.off()

#+
number_of_pvalues_formatted <- prettyNum(nrow(number_of_pvalues), big.mark=',')
median_number_of_pvalues_formatted <- prettyNum(round(median(number_of_pvalues$n_pvalues)), big.mark=',')

fig_cap <- glue("__Reduced number of features in anti-conservative and uniform p value sets.__ 
                (__A__) P value set size distribution. Dashed line denotes the median number of features ({median_number_of_pvalues_formatted}). From each GEO series only one random set was considered, N = {number_of_pvalues_formatted} p value sets. 
                (__B__) Robust linear modeling of number of features in anti-conservative and uniform vs. non-anti-conservative p value sets [log10_n_pvalues ~ anticons, student’s t likelihood], N = {number_of_pvalues_formatted}. 
                Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible region, respectively.
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/log10_n_pvalues%20~%20anticons.rds."
                )

#+ fig.cap=fig_cap
p


#+ s2fig
f <- anticons ~ year
data <- pvalues_sample
family <- bernoulli()
conditions <- NULL
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_year.rds"),
           file_refit = "on_change")
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions),
          line_args = list(color = "black"),
          plot = FALSE)

fig_cap <- glue("__The increasing proportion of anti-conservative histograms.__ 
                Bernoulli model [{deparse(f)}], N = {prettyNum(summary(mod)$nobs, big.mark=',')}. 
                Lines denote best fit of linear model. Shaded area denotes 95% credible region.
                The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_year.rds.")

#+ fig.cap=fig_cap
py <- p$year + 
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Year")
py
#+ include=FALSE
# tiff(here("figures/S2Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# py
# dev.off()

#' 
#+ s3fig
# f <- anticons ~ year + (year | de_tool)
# mod <- brm(formula = f, 
#            data = data, 
#            family = family, 
#            chains = chains, 
#            cores = cores, 
#            refresh = refresh,
#            iter = ifelse(is_ci(), 400, 2000),
#            control = list(adapt_delta = 0.99, max_treedepth = 12),
#            file = here("results/models/anticons_year__year_detool.rds"),
#            file_refit = "on_change"
# )
# 
# draws <- data %>% 
#   add_epred_draws(mod, ndraws = 100)
# 
# p <- draws %>% 
#   ggplot(aes(year, .epred, group = 1)) +
#   stat_lineribbon(point_interval = median_qi, .width = 0.95, fill = "gray65")
# 
# fig_cap <- glue("__Trends of anti-conservative p value histograms across differential expression analysis tools__, 
#                 two-level binomial logistic model [{deparse(f)}], 
#                 N = {prettyNum(summary(mod)$nobs, big.mark=',')}.
#                 Lines denote best fit of linear model. Shaded area denotes 95% credible region.
#                 The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_year__year_detool.rds.")
# 
# #+fig.cap=fig_cap
# pat <- p + 
#   scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
#   labs(y = "Proportion of anti-conservative\np value histograms",
#        x = "Year") +
#   facet_wrap(~ de_tool, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# pat
data <- pvalues_sample %>% 
  mutate_at("year", ~.x - min(.x)) %>%
  select(Class, year, de_tool) %>% 
  drop_na()

priors <- c(
  set_prior("lkj(3)", class = "cor"),
  set_prior("student_t(3, 0, 1)", class = "Intercept"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "muuniform"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "muuniform"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "mubimodal"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "mubimodal"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "muconservative"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "muconservative"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "muother"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "muother")
)

fa <- Class ~ year + (year | de_tool)
moda <- brm(
  formula = fa, 
  data = data, 
  family = categorical(), 
  chains = chains, 
  cores = chains, 
  refresh = refresh,
  prior = priors,
  iter = ifelse(is_ci(), 400, 2400),
  file = here("results/models/Class_year__year_detool_year.rds"),
  file_refit = "never"
)

draws <- moda$data %>% 
  select(year, de_tool) %>% 
  filter(de_tool != "seurat") %>% 
  distinct() %>% 
  add_epred_draws(moda)
p2a <- draws %>% 
  ggplot(aes(year + min(pvalues_sample$year), .epred)) +
  stat_summary(fun.data = mean_hdci, aes(fill = .category), geom = "ribbon",  alpha = 0.2, fun.args = list(.width = 0.95)) +
  stat_summary(fun = mean, aes(color = .category), geom = "line") +
  facet_wrap(~str_replace(de_tool, "-", "-\n"), nrow = 2, scales = "free_x") +
  labs(y = "Proportion",
       x = "Year") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_colour_manual(values = c("#00BF7D", "#A3A500", "#F8766D", "#00B0F6", "#E76BF3")) +
  scale_fill_manual(values = c("#00BF7D", "#A3A500", "#F8766D", "#00B0F6", "#E76BF3"))

#'
#'
#+
data <- pvalues_sample %>% 
  filter(year >= 2018) %>% # keep only values from 2018 and up!
  count(Class, de_tool) %>% 
  complete(Class, de_tool, fill = list(n = 0)) %>% 
  group_by(de_tool) %>% 
  mutate(total_in_de_tool = sum(n)) %>% 
  ungroup() %>% 
  arrange(de_tool)

priors <- c(
  set_prior("student_t(3, 0, 25)", class = "b")
)

#+
fb <- n | trials(total_in_de_tool) ~ Class + de_tool + Class:de_tool
modb <- brm(
  formula = fb,
  data = data,
  family = binomial(),
  prior = priors,
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  control = list(max_treedepth = 15),
  iter = ifelse(is_ci(), 400, 4000),
  file = here("results/models/n__trials(total_in_de_tool)__Class_de_tool_Class:de_tool_2018up.rds"),
  file_refit = "on_change"
)

draws <- modb$data %>% 
  data_grid(Class, de_tool, total_in_de_tool = 1000) %>% 
  add_linpred_draws(modb) %>% 
  mutate_at(".linpred", inv_logit_scaled)
p2b <- draws %>% 
  ungroup() %>% 
  filter(Class != "uniform") %>% 
  ggplot(aes(de_tool, .linpred)) +
  stat_pointinterval(point_size = 1) +
  facet_wrap(~Class, nrow = 1) +
  labs(y = "Proportion") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )

p2_full <- p2a / p2b + plot_annotation(tag_levels = "A") +  plot_layout(guides = 'auto', heights = c(2, 1))

fig_cap <- glue(
  "__Association of the p-value histogram class with a differential expression analysis tool.__ 
(__A__) Time courses for proportions of different p-value histogram classes for the nine most frequent DE analysis platforms. 
Lines denote best fit of the model [{deparse(fa)}, categorical likelihood]. 
Shaded areas denote 95% credible regions. N = {prettyNum(summary(moda)$nobs, big.mark=',')}. 
(__B__) Association of p-value histogram type with DE analysis tool; data is restricted to 2018-2020 GEO submissions. 
Points denote best fit of the model [{deparse(fb)}, binomial likelihood]. 
Thick and thin lines denote 66% and 95% credible intervals, respectively. N = {prettyNum(sum(data$n), big.mark=',')}. 
The model object related to panel A can be downloaded from 
https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/Class_year__year_detool_year.rds. 
The model object related to panel B can be downloaded from 
https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/n__trials(total_in_de_tool)__Class_de_tool_Class:de_tool_2018up.rds."
)

#+fig.cap=fig_cap, fig.height=8
p2_full

#+ include=FALSE
# tiff(here("figures/S3Fig.tiff"), height = 12, width = 18, units = "cm", res = 300)
# p2_full
# dev.off()

#' 
#+ s4fig
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  select(Accession, id, Set, anticons, year, model) %>% 
  drop_na()
f <- anticons ~ year + (year | model)
mod <- brm(formula = f, 
           data = data, 
           family = family,
           iter = ifelse(is_ci(), 400, 4000),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_year__year_model.rds"),
           file_refit = "on_change")

draws <- data %>% 
  filter(
    model %in% c("illumina hiseq 4000",          "nextseq 500",                  "illumina miseq",               "illumina hiseq 2500",         
                      "illumina hiseq 2000",          "illumina hiseq 1000",          "illumina genome analyzer iix", "ion torrent proton",          
                      "illumina hiseq 1500",          "ab 5500 genetic analyzer",     "illumina hiseq 3000",          "bgiseq-500",                  
                      "hiseq x ten",                  "illumina hiscansq",            "illumina genome analyzer",     "ab solid 4 system",           
                      "nextseq 550",                  "illumina genome analyzer ii",  "ab 5500xl genetic analyzer",   "illumina novaseq 6000",       
                      "ion torrent pgm",              "illumina hiseq x ten",         "ab solid system 3.0",          "ab solid system")
    ) %>% 
  select(model, year) %>%
  distinct() %>% 
  add_epred_draws(mod)

p <- draws %>% 
  ggplot(aes(year, .epred, group = 1)) +
  stat_lineribbon(point_interval = median_qi, .width = 0.95, fill = "gray65") +
  facet_wrap(~ model, labeller = label_wrap_gen(width = 18))

fig_cap <- glue("__All sequencing instrument models are associated with temporally increasing anti-conservative p value histograms__,
                two-level binomial logistic model [{deparse(f)}], N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Only GEO submissions utilizing 
                single sequencing platform were used for model fitting. Lines denote best fit of linear model. Shaded area denotes 95% credible region. 
                The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_year__year_model.rds.")

#+ fig.height=8, fig.cap=fig_cap, include=FALSE
pim <- p + 
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(y = "Proportion of anti-conservative p value histograms",
       x = "Year") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 8))
pim

#+ include=FALSE
# tiff(here("figures/S4Fig.tiff"), height = 23, width = 18, units = "cm", res = 300)
# pim
# dev.off()

#'
#+ s5fig
p <- pvalues_sample %>% 
  count(year, de_tool) %>% 
  add_count(year, name = "total", wt = n) %>% 
  mutate(Proportion = n / total) %>% 
  ggplot() +
  geom_line(aes(year, Proportion, color = de_tool), size = 1) +
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(x = "Year") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3))

fig_cap <- glue("__No single differential expression analysis tool dominates the field.__ 
                Y-axis shows the proportion of analysis platforms, 
                x-axis shows publication year of GEO submission, N = {prettyNum(sum(p$data$n), big.mark=',')}.")

#+ fig.cap=fig_cap
p
#+ include=FALSE
# tiff(here("figures/S5Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#+ s6figa
y_title <- "Prop. anti-cons."
f <- anticons ~ de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_detool.rds"),
           file_refit = "on_change")

draws_6a <- data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod)
pa <- draws_6a %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

#'
#+ s6figb
data <- pvalues
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_detool_all.rds"),
           file_refit = "on_change")

pb <- data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)

#' 
#+ s6figc
f <- anticons ~ year + de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_year_detool.rds"),
           file_refit = "on_change")

pc <- data %>% 
  data_grid(de_tool, year = 2020) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pcN <- prettyNum(summary(mod)$nobs, big.mark=',')
pcF <- deparse(f)

#' 
#+ s6figd
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == 9606 ~ "human",
    tax_id == 10090 ~ "mouse",
    TRUE ~ "other"
  ))
f <- anticons ~ organism + de_tool
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_organism_detool.rds"),
           file_refit = "on_change")

pd <- data %>% 
  data_grid(de_tool, organism = "human") %>%
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdN <- prettyNum(summary(mod)$nobs, big.mark=',')
pdF <- deparse(f)

#' 
#+ s6fige
f <- anticons ~ de_tool + (1 | model)
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  drop_na(de_tool, model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_detool__1_model.rds"),
           file_refit = "on_change")

draws <- as.data.frame(mod)
pe <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(data$de_tool))) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

#' 
#+ s6figf
f <- anticons ~ de_tool + (de_tool | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_detool__detool_model.rds"),
           file_refit = "on_change")

draws <- as.data.frame(mod)
pf <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(data$de_tool))) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pfN <- prettyNum(summary(mod)$nobs, big.mark=',')
pfF <- deparse(f)

fig_cap <- glue('__DE analysis tool conditional effects from binomial logistic models for proportion of anti-conservative p value histograms.__ 
                (__A__) Simple model [{paF}], N = {paN}. 
                (__B__) Simple model [{pbF}] fitted on complete data, N = {pbN}.
                (__C__) Model conditioned on year of GEO submission [{pcF}], N = {pcN}.
                (__D__) Model conditioned on studied organism (human/mouse/other) [{pdF}], N = {pdN}. 
                (__E__) Varying intercept model [{peF}] where "model" stands for sequencing instrument model, N = {peN}. 
                (__F__) Varying intercept and slope model [{pfF}], N = {pfN}. 
                Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. 
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool_all.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_year_detool.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_organism_detool.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool__1_model.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool__detool_model.rds.')

#+ s6fig, fig.cap=fig_cap
pde <- (pa + pb + pc) / (pd + pe + pf) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 10, hjust = 1, vjust = -1), plot.margin = unit(c(0.5,0.2,0.2,0.2), "cm"))
pde

#+ include=FALSE
# tiff(here("figures/S6Fig.tiff"), height = 14, width = 18, units = "cm", res = 300)
# pde
# dev.off()


#'
#'
#+
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == 9606 ~ "human",
    tax_id == 10090 ~ "mouse",
    TRUE ~ "other"
  )) %>% 
  inner_join(n_samples) %>% 
  drop_na(organism, N, de_tool)
f <- anticons ~ organism + (N | de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = bernoulli(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons ~ organism + (N | de_tool).rds"),
           file_refit = "on_change")


f <- pi0 ~ organism + (N | de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = Beta(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0 ~ organism + (N | de_tool).rds"),
           file_refit = "on_change")

#' 
#+ s7figa
f <- pi0 ~ de_tool
family <- Beta()

#####
mod <- read_rds(here("results/models/pi0_detool_sample.rds"))
draws_7aa <- mod$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(mod)
paa <- draws_7aa %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

paaN <- prettyNum(summary(mod)$nobs, big.mark=',')
paaF <- deparse(f)

#####

data <- pvalues %>% 
  drop_na(pi0, de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0_detool_full_data.rds"),
           file_refit = "on_change")

draws_7a <- data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod)
pa <- draws_7a %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

#' 
#+ s7figb
f <- pi0 ~ year + de_tool
data <- pvalues_sample %>% 
  drop_na(pi0, year, de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0_year_detool.rds"),
           file_refit = "on_change")

pb <- data %>% 
  data_grid(de_tool, year = 2020) %>%
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)

#' 
#+ s7figc
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == 9606 ~ "human",
    tax_id == 10090 ~ "mouse",
    TRUE ~ "other"
  )) %>% 
  drop_na(pi0, organism, de_tool)
f <- pi0 ~ organism + de_tool
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("results/models/pi0_organism_detool.rds"),
           file_refit = "on_change")

pc <- data %>% 
  data_grid(de_tool, organism = "human") %>%
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pcN <- prettyNum(summary(mod)$nobs, big.mark=',')
pcF <- deparse(f)

#' 
#+ s7figd
f <- pi0 ~ de_tool + (1 | model)
mod <- brm(formula = f, 
           data = data %>% drop_na(model), 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0_detool__1_model.rds"),
           file_refit = "on_change"
           )

draws <- as.data.frame(mod)
pd <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(data$de_tool))) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdN <- prettyNum(summary(mod)$nobs, big.mark=',')
pdF <- deparse(f)

#' 
#+ s7fige
f <- pi0 ~ de_tool + (de_tool | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0_detool__detool_model.rds"),
           file_refit = "on_change")

draws <- as.data.frame(mod)
pe <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(data$de_tool))) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

p <- (paa + pa + pb) / (pc + pd + pe) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))
fig_cap <- glue("__DE analysis tool conditional effects from beta regression modeling of $\\pi_0$.__ 
                (__A__) Simple model [{paaF}] fitted on sample, N = {paaN}. 
                (__B__) Simple model [{paF}] fitted on complete data, N = {paN}. 
                (__C__) Model conditioned on year of GEO submission [{pbF}], N = {pbN}. 
                (__D__) Model conditioned on studied organism (human/mouse/other) [{pcF}], N = {pcN}. 
                (__E__) Varying intercept model [{pdF}] where 'model' stands for sequencing instrument model, N = {pdN}. 
                (__F__) Varying intercept/slope model [{peF}], N = {peN}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively.
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool_sample.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool_full_data.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_year_detool.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_organism_detool.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool__1_model.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool__detool_model.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S7Fig.tiff"), height = 14, width = 18, units = "cm", res = 300)
# p
# dev.off()

#'
#' 
#+
parsed_suppfiles_rerun <- read_csv(here("results/data/parsed_suppfiles_rerun.csv")) %>% 
  distinct() %>% 
  filter(Type == "raw") %>% 
  select(-Type)
pvalues_rerun <- pvalues %>% 
  inner_join(parsed_suppfiles_rerun, by = c("id", "Set"), suffix = c(".orig", ".rerun")) %>% 
  rename_all(str_remove, ".orig") %>% 
  rename(prop_alpha_pval = prop_FDR_pval, alpha_pval = FDR_pval)

pi0methodsa <- pvalues_rerun %>% 
  filter(!is.na(pi0)) %>% 
  select(Accession, id, pi0, pi0.rerun) %>% 
  pivot_longer(cols = starts_with("pi0")) %>% 
  mutate(name = ifelse(name == "pi0", "Local FDR", "Smoother")) %>% 
  ggplot() +
  geom_density(aes(value, color = name)) +
  labs(x = expression(pi[0]), y = "Density") +
  theme(legend.title = element_blank())

pi0methodsb <- pvalues_rerun %>% 
  filter(!is.na(pi0)) %>% 
  ggplot() +
  geom_point(aes(pi0, pi0.rerun), size = 1/20) +
  geom_abline(linetype = "dashed") +
  labs(x = "Local FDR", y = "Smoother")

pi0methodsplot <- pi0methodsa + pi0methodsb + plot_layout(widths = c(3, 2)) + plot_annotation(tag_levels = "A")
fig_cap <- glue("__Comparison of $\\pi_0$ values computed by two different methods.__ Local FDR method is from limma R package [@limma] function propTrueNull. Smoother method is from qvalue R package [@qvalue]. A, density histogram. B, scatter plot. Dashed line has intercept = 0 and slope = 1.")

#+ fig.cap=fig_cap
pi0methodsplot

#+ include=FALSE
# tiff(here("figures/S8Fig.tiff"), height = 9, width = 18, units = "cm", res = 300)
# pi0methodsplot
# dev.off()

#' 
#+ s8fig
f <- pi0 ~ model
data <- data %>% 
  drop_na(pi0, model)
mod <- brm(formula = f, 
           data = data, 
           family = Beta(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__model.rds"),
           file_refit = "on_change")
p <- data %>% 
  data_grid(model) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, model)) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0])) +
  theme(axis.title.y = element_blank()) +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on sequencing instrument model.__ Points denote best fit of linear model ([{deparse(f)}], beta distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0__model.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S9Fig.tiff"), height = 12, width = 18, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s9fig
f <- pi0 ~ library_strategy
mod <- brm(formula = f, 
           data = data, 
           family = Beta(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__librarystrategy.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(library_strategy) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, library_strategy)) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0]), y = "Library strategy") +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on library strategy.__ Points denote best fit of linear model ([{deparse(f)}], beta distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0__librarystrategy.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S10Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s10fig
f <- pi0 ~ library_selection
mod <- brm(formula = f, 
           data = data, 
           family = Beta(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__libraryselection.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(library_selection) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, library_selection)) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0]), y = "Library selection") +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on library selection.__ Points denote best fit of linear model ([{deparse(f)}, beta likelihood], N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0__libraryselection.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S11Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s11fig
f <- pi0 ~ library_layout
mod <- brm(formula = f, 
           data = data, 
           family = Beta(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores,
           iter = ifelse(is_ci(), 400, 4000),
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__1_librarylayout.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(library_layout) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, library_layout)) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0]), y = "Library layout") +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on library layout.__ Points denote best fit of linear model ([{deparse(f)}, beta likelihood], N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0__librarylayout.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S12Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s12fig
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  drop_na(model)
f <- anticons ~ model
mod <- brm(formula = f, 
           data = data, 
           family = bernoulli(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__1_model.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(model) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, model)) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms") +
  theme(axis.title.y = element_blank()) +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on sequencing platform.__ Points denote best fit of linear model ([{deparse(f)}, bernoulli likelihood], N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons__model.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S13Fig.tiff"), height = 12, width = 18, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s13fig
f <- anticons ~ library_strategy
mod <- brm(formula = f, 
           data = data, 
           family = bernoulli(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__librarystrategy.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(library_strategy) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, library_strategy)) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms", y = "Library strategy") +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on library strategy.__ Points denote best fit of linear model ([{deparse(f)}, bernoulli likelihood], N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons__librarystrategy.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S14Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s14fig
f <- anticons ~ library_selection
mod <- brm(formula = f, 
           data = data, 
           family = bernoulli(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__libraryselection.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(library_selection) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, library_selection)) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms", y = "Library selection") +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on library selection.__ Points denote best fit of linear model ([{deparse(f)}, bernoulli likelihood], N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons__libraryselection.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S15Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#' 
#+ s15fig
f <- anticons ~ library_layout
mod <- brm(formula = f, 
           data = data, 
           family = bernoulli(), 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__librarylayout.rds"),
           file_refit = "on_change")

p <- data %>% 
  data_grid(library_layout) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, library_layout)) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms", y = "Library layout") +
  scale_x_continuous(limits = c(0, 1))

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on library layout.__ Points denote best fit of linear model ([{deparse(f)}, bernoulli likelihood], N = {prettyNum(summary(mod)$nobs, big.mark=',')}). Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons__librarylayout.rds.")

#+ fig.cap=fig_cap
p

#+ include=FALSE
# tiff(here("figures/S16Fig.tiff"), height = 7, width = 10, units = "cm", res = 300)
# p
# dev.off()

#'
#'
#+
pvalues_filtered_sample <- pvalues_filtered %>% 
  inner_join(pvalues_sample %>% select(Accession, id, Set)) %>% 
  left_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == "9606" ~ "human",
    tax_id == "10090" ~ "mouse",
    TRUE ~ "other"
  ))

#+ s16figa
y_title <- "Prop. anti-cons."
f <- anticons ~ de_tool
data <- pvalues_filtered_sample
mod <- brm(
  formula = f, 
  data = data, 
  family = bernoulli(), 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/anticons_detool_filtered.rds"),
  file_refit = "on_change"
)

draws_16a <- data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod)
pa <- draws_16a %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

es_anticons_prop <- draws_6a %>% 
  ungroup() %>% 
  rename(raw = .epred) %>% 
  select(de_tool, .draw, raw) %>% 
  left_join(
    draws_16a%>% 
      ungroup() %>% 
      rename(filtered = .epred) %>% 
      select(de_tool, .draw, filtered)) %>% 
  mutate(es = filtered - raw)

pg <- es_anticons_prop %>% 
  ggplot(aes(es, de_tool)) +
  stat_halfeye(point_size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size") +
  scale_y_discrete(limits = rev) +
  theme(axis.title.y = element_blank())

#'
#+ s16figb
mod <- brm(
  formula = f, 
  data = pvalues_filtered, 
  family = bernoulli(), 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/anticons_detool_all_filtered.rds"),
  file_refit = "on_change"
)

pb <- data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)

#' 
#+ s16figc
mod <- brm(
  formula = anticons ~ year + de_tool, 
  data = pvalues_filtered_sample, 
  family = bernoulli(), 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/anticons_year_detool_filtered.rds"),
  file_refit = "on_change"
)

pc <- data %>% 
  data_grid(de_tool, year = 2020) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pcN <- prettyNum(summary(mod)$nobs, big.mark=',')
pcF <- deparse(f)

#' 
#+ s16figd
f <- anticons ~ organism + de_tool
mod <- brm(
  formula = f, 
  data = pvalues_filtered_sample, 
  family = bernoulli(), 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/anticons_organism_detool_filtered.rds"),
  file_refit = "on_change"
)

pd <- data %>% 
  data_grid(de_tool, organism = "human") %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  # scale_y_continuous(limits = c(0, 1)) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdN <- prettyNum(summary(mod)$nobs, big.mark=',')
pdF <- deparse(f)

#' 
#+ s16fige
f <- anticons ~ de_tool + (1 | model)
mod <- brm(
  formula = f, 
  data = data, 
  family = bernoulli(), 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = here("results/models/anticons_detool__1_model_filtered.rds"),
  file_refit = "on_change"
)

draws <- as.data.frame(mod)
pe <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(data$de_tool))) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

#' 
#+ s16figf
f <- anticons ~ de_tool + (de_tool | model)
mod <- brm(
  formula = f, 
  data = data, 
  family = bernoulli(), 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = here("results/models/anticons_detool__detool_model_filtered.rds"),
  file_refit = "on_change"
)

draws <- as.data.frame(mod)
pf <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(data$de_tool))) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pfN <- prettyNum(summary(mod)$nobs, big.mark=',')
pfF <- deparse(f)

fig_cap <- glue('__DE analysis tool conditional effects from binomial logistic models for proportion of anti-conservative p value histograms after filtering to remove p values for lowly expressed features.__ 
                Complete data contains all observations where expression level data was available. Sample contains observations from the original unfiltered sample where filtering was possible.
                (__A__) Simple model [{paF}], N = {paN}. 
                (__B__) Simple model [{pbF}] fitted on complete data, N = {pbN}. 
                (__C__) Model conditioned on year of GEO submission [{pcF}], N = {pcN}. 
                (__D__) Model conditioned on studied organism (human/mouse/other) [{pdF}], N = {pdN}. 
                (__E__) Varying intercept model [{peF}] where "model" stands for sequencing instrument model, N = {peN}. 
                (__F__) Varying intercept and slope model [{pfF}], N = {pfN}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively.
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool_filtered.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool_all_filtered.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_year_detool_filtered.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_organism_detool_filtered.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool__1_model_filtered.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/anticons_detool__detool_model_filtered.rds.')

#+ s16fig
p16 <- (pa + pb + pc) / (pd + pe + pf) +  
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

#+ include=FALSE
# tiff(here("figures/S17Fig.tiff"), height = 15, width = 18, units = "cm", res = 300)
# p16
# dev.off()

#' 
#+ s17figaa
f <- pi0 ~ de_tool
family <- Beta()
mod <- brm(
  formula = f, 
  data = pvalues_filtered_sample, 
  family = family, 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/pi0_detool_sample_filtered.rds"),
  file_refit = "on_change"
)

draws_17aa <- mod$data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod)
paa <- draws_17aa %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

paaN <- prettyNum(summary(mod)$nobs, big.mark=',')
paaF <- deparse(f)

pi0_detool_sample <- read_rds(here("results/models/pi0_detool_sample.rds"))
es_pi0_prop <- pi0_detool_sample$data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(pi0_detool_sample) %>% 
  ungroup() %>% 
  rename(raw = .epred) %>% 
  select(de_tool, .draw, raw) %>% 
  left_join(
    draws_17aa %>% 
      ungroup() %>% 
      rename(filtered = .epred) %>% 
      select(de_tool, .draw, filtered)) %>% 
  mutate(es = filtered - raw)

pg <- es_pi0_prop %>% 
  ggplot(aes(es, de_tool)) +
  stat_halfeye(point_size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size") +
  scale_y_discrete(limits = rev) +
  theme(axis.title.y = element_blank())

#' 
#+ s17figa
f <- pi0 ~ de_tool
family <- Beta()
mod <- brm(
  formula = f, 
  data = pvalues_filtered, 
  family = family, 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/pi0_detool_full_data_filtered.rds"),
  file_refit = "on_change"
)

draws_17a <- mod$data %>% 
  data_grid(de_tool) %>% 
  add_epred_draws(mod)
pa <- draws_17a %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

#' 
#+ s17figb
f <- pi0 ~ year + de_tool
mod <- brm(
  formula = f, 
  data = pvalues_filtered_sample, 
  family = family, 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/pi0_year_detool_filtered.rds"),
  file_refit = "on_change"
)

pb <- mod$data %>% 
  data_grid(de_tool, year = 2020) %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)

#' 
#+ s17figc
f <- pi0 ~ organism + de_tool
mod <- brm(
  formula = f, 
  data = pvalues_filtered_sample, 
  family = family, 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  file = here("results/models/pi0_organism_detool_filtered.rds"),
  file_refit = "on_change"
)

pc <- mod$data %>% 
  data_grid(de_tool, organism = "human") %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pcN <- prettyNum(summary(mod)$nobs, big.mark=',')
pcF <- deparse(f)

#' 
#+ s17figd
f <- pi0 ~ de_tool + (1 | model)
mod <- brm(
  formula = f, 
  data = pvalues_filtered_sample, 
  family = family, 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = here("results/models/pi0_detool__1_model_filtered.rds"),
  file_refit = "on_change"
)

draws <- as.data.frame(mod)
pd <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(pvalues_filtered_sample$de_tool))) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdN <- prettyNum(summary(mod)$nobs, big.mark=',')
pdF <- deparse(f)

#' 
#+ s17fige
f <- pi0 ~ de_tool + (de_tool | model)
mod <- brm(
  formula = f, 
  data = pvalues_filtered_sample, 
  family = family, 
  prior = prior("student_t(3, 0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = here("results/models/pi0_detool__detool_model_filtered.rds"),
  file_refit = "on_change"
)

draws <- as.data.frame(mod)
pe <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  scale_x_discrete(labels = sort(unique(pvalues_filtered_sample$de_tool))) +
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

p17 <- (paa + pa + pb) / (pc + pd + pe) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

fig_cap <- glue("__DE analysis tool conditional effects from beta regression models of $\\pi_0$ after filtering to remove p values for lowly expressed features.__ 
Complete data contains all observations where expression level data was available. Sample contains observations from the original unfiltered sample where filtering was possible.
                (__A__) Simple model [{paaF}], N = {paaN}. 
                (__B__) Simple model [{paF}] fitted on complete data, N = {paN}. 
                (__C__) Model conditioned on year of GEO submission [{pbF}], N = {pbN}. 
                (__D__) Model conditioned on studied organism (human/mouse/other) [{pcF}], N = {pcN}. 
                (__E__) Varying intercept model [{pdF}] where 'model' stands for sequencing instrument model, N = {pdN}. 
                (__F__) Varying intercept/slope model [{peF}], N = {peN}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively.
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool_sample_filtered.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool_full_data_filtered.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_year_detool_filtered.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_organism_detool_filtered.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool__1_model_filtered.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/src/v0.1/models/pi0_detool__detool_model_filtered.rds.")

#+ include=FALSE
# tiff(here("figures/S18Fig.tiff"), height = 15, width = 18, units = "cm", res = 300)
# p17
# dev.off()

#'
#'
#+
anticons_det <- read_rds(here("results/models/anticons_detool.rds"))
anticons_det_filt <- read_rds(here("results/models/anticons_detool_filtered.rds"))
pi0_det <- read_rds(here("results/models/pi0_detool_sample.rds"))
pi0_det_filt <- read_rds(here("results/models/pi0_detool_sample_filtered.rds"))

de_tools <- pvalues_sample %>% 
  select(de_tool) %>% 
  distinct()

draws_anticons_det <- de_tools %>% 
  add_epred_draws(anticons_det, value = "raw") %>% 
  ungroup() %>% 
  select(de_tool, .draw, raw)
draws_anticons_det_filt <- anticons_det_filt$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(anticons_det_filt, value = "filtered") %>% 
  ungroup() %>% 
  select(de_tool, .draw, filtered)
draws_pi0_det <- pi0_det$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(pi0_det, value = "raw") %>% 
  ungroup() %>% 
  select(de_tool, .draw, raw)
draws_pi0_det_filt <- pi0_det_filt$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(pi0_det_filt, value = "filtered") %>% 
  ungroup() %>% 
  select(de_tool, .draw, filtered)

draws_ac_merged <- draws_anticons_det %>% 
  inner_join(draws_anticons_det_filt)

pd <- position_dodge(0.6)
pa <- draws_ac_merged %>% 
  pivot_longer(cols = c("raw", "filtered")) %>% 
  mutate(name = factor(name, levels = c("raw", "filtered"))) %>% 
  ggplot(aes(value, de_tool)) +
  stat_pointinterval(aes(color = name), point_size = 1, position = pd) +
  labs(x = "Prop. anti-cons.") +
  scale_y_discrete(limits = rev) +
  scale_color_discrete("P value\nset") +
  theme(axis.title.y = element_blank())

pb <- draws_ac_merged %>% 
  mutate(es = filtered - raw) %>% 
  ggplot(aes(es, de_tool)) +
  stat_pointinterval(point_size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1/3) +
  labs(x = "Effect size") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-0.1, 0.45)) +
  theme(axis.title.y = element_blank())

draws_pi0_merged <- draws_pi0_det %>% 
  inner_join(draws_pi0_det_filt)

pc <- draws_pi0_merged %>% 
  pivot_longer(cols = c("raw", "filtered")) %>% 
  mutate(name = factor(name, levels = c("raw", "filtered"))) %>% 
  ggplot(aes(value, de_tool)) +
  stat_pointinterval(aes(color = name), point_size = 1, position = pd) +
  labs(x = expression(pi[0])) +
  scale_y_discrete(limits = rev) +
  scale_color_discrete("P value\nset") +
  theme(axis.title.y = element_blank())

pd <- draws_pi0_merged %>% 
  mutate(es = filtered - raw) %>% 
  ggplot(aes(es, de_tool)) +
  stat_pointinterval(point_size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1/3) +
  labs(x = "Effect size") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-0.1, 0.45)) +
  theme(axis.title.y = element_blank())

patchwork <- (pa + pb) / (pc + pd) + plot_layout(nrow = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A")

paN <- prettyNum(summary(anticons_det_filt)$nobs, big.mark=',')
paF <- deparse(anticons ~ de_tool)
pcN <- prettyNum(summary(pi0_det_filt)$nobs, big.mark=',')
pcF <- deparse(pi0 ~ de_tool)


fig_cap <- glue(
  "__Removal of low-count features results in an increasing proportion of anti-conservative p-value histograms.__
  (__A__) Anti-conservative p-value histogram proportions in raw and filtered p-value sets for DE analysis programs. 
  Raw p-value data is the same as in S5 Fig. A. Filtered p-value data is from a simple Bernoulli model [{paF}], N = {paN}. 
  (__B__) Effect size of low-count feature filtering to proportion of anti-conservative p-values. 
  (__C__) $\\pi_0$ estimates for raw and filtered p-value sets. 
  Raw p-value data is the same as in S6 Fig. A and filtered p-value data is from the beta model [{pcF}], N = {pcN}. 
  (__D__) Effect size of low-count feature filtering to $\\pi_0$.
  Points denote model best fit. Thick- and thin lines denote 66% and 95% CIs, respectively."
  )

#+ fig.cap=fig_cap
patchwork

#####

#' 
#+ s18fig
qc_threshold <- function(n_pvals, bins, fdr) {
  qbinom(1 - 1 / bins * fdr, n_pvals, 1 / bins)
}

qc_thr <- qc_threshold(n_pvals = 20000, bins = 40, fdr = 0.05)

p <- de_simulation_results %>% 
  mutate(
    reps = factor(reps, levels = c("N=3", "N=6", "N=10")),
    n_eff = factor(n_eff, levels = c("Effects=100", "Effects=200", "Effects=400", "Effects=800"))
  ) %>% 
  ggplot() +
  geom_histogram(aes(pvalue), binwidth = 1/40, center = 1/80) +
  geom_hline(yintercept = qc_thr, color = "red") +
  facet_grid(reps~n_eff) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

#+ include=FALSE
# tiff(here("figures/S19Fig.tiff"), height = 8, width = 12, units = "cm", res = 300)
# p
# dev.off()

#+
fig_cap <- glue('__Simulated RNA-seq data shows that histograms from p value sets with around one hundred 
                 true effects out of 20,000 features can be classified as "uniform".__
                RNA-seq data was simulated with polyester R package [@frazee2015] on 20,000 transcripts from human transcriptome 
                 using grid of 3, 6, and 10 replicates and 100, 200, 400, and 800 effects for two groups. 
                 Fold changes were set to 0.5 and 2.
                Differential expression was assessed using DESeq2 R package [@love2014] using default settings 
                and group 1 versus group 2 contrast. 
                Effects denotes in facet labels the number of true effects and N denotes number of replicates.
                Red line denotes QC threshold used for dividing p histograms into discrete classes.
                Code and workflow used to run these simulations is available on Github: https://github.com/rstats-tartu/simulate-rnaseq.
                Raw data of the figure is available on Zenodo https://zenodo.org with doi: 10.5281/zenodo.4463803.')

#+ fig.cap=fig_cap
p


#' \newpage
#'
#'
#' ## References
#'


