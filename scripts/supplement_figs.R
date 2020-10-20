#' ---
#' title: Supplementary figures
#' date: ""
#' author: ""
#' output:
#'    bookdown::pdf_book:
#'        number_sections: FALSE
#'        toc: FALSE
#' ---

#+ include=FALSE
knitr::opts_chunk$set(echo = TRUE, message = FALSE, comment = FALSE, warning = FALSE)


#+ libs
library(readr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(stringr)
library(cowplot)
library(patchwork)
library(brms)
library(tidybayes)
library(here)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

#'
#+ params
chains <- 4
cores <- chains
refresh = 0

#+ data
pvalues <- read_csv(here("output/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample <- read_csv(here("output/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
sequencing_metadata <- read_csv(here("output/sequencing_metadata_unique_platform.csv"))

#+ fig.cap = "The increasing proportion of anti-conservative histograms. Binomial logistic model: $anticons \\sim year$, N = 2109."
f <- anticons ~ year
family <- bernoulli()
data <- pvalues_sample
conditions <- NULL
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/anticons_year.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions),
          plot = FALSE)
p$year + 
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(y = "Proportion of anti-conservative p histograms",
       x = "Year")

#' 
#+ fig.cap="A 2-level binomial logistic model $anticons \\sim year + (year | DEtool)$ reveals that all differential expression analysis tools are associated with temporally increasing anti-conservative p histograms, N = 2109."
f <- anticons ~ year + (year | de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons_year__year_detool.rds"))
conditions <- make_conditions(data, vars = "de_tool")
row.names(conditions) <- conditions$de_tool
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions, 
                              re_formula = NULL),
          plot = FALSE)
p$year + 
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(y = "Proportion of anti-conservative p histograms",
       x = "Year")

#' 
#+ FigS3, fig.cap="A 2-level binomial logistic model $anticons \\sim year + (year  | model)$ reveals that all sequencing instrument models are associated with temporally increasing anti-conservative p histograms, N = 1718. Only GEO submissions utilizing single sequencing platform were used for model fitting."
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata)
f <- anticons ~ year + (year | model)
mod <- brm(formula = f, 
           data = data, 
           family = family,
           iter = 4000,  
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons_year__year_model.rds"))
conditions <- make_conditions(data, vars = "model")
row.names(conditions) <- str_wrap(conditions$model, width = 20)
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions, 
                              re_formula = NULL),
          plot = FALSE)
p$year + 
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(y = "Proportion of anti-conservative p histograms",
       x = "Year")

#'
#+ FigS4, fig.cap="No single data analysis platform dominates the field. Y-axis shows the proportion of analysis platforms, x-axis shows publication year of GEO submission, N = 1733."
data %>% 
  count(year, de_tool) %>% 
  add_count(year, name = "total", wt = n) %>% 
  mutate(Proportion = n / total) %>% 
  ggplot() +
  geom_line(aes(year, Proportion, color = de_tool), size = 1) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(x = "Year") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")


#+ FigS5a
f <- anticons ~ de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/anticons_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pa <- p$de_tool + 
  labs(y = "Proportion of anti-conservative\np histograms",
       x = "DE analysis tool")

#'
#+ FigS5b
data <- pvalues
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/anticons_detool_all.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pb <- p$de_tool + 
  labs(y = "Proportion of anti-conservative p histograms",
       x = "DE analysis tool") + 
  theme(axis.title.y = element_blank())

#' 
#+ FigS5c
f <- anticons ~ year + de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/anticons_year_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pc <- p$de_tool + 
  labs(y = "Proportion of anti-conservative p histograms",
       x = "DE analysis tool") + 
  theme(axis.title.y = element_blank())

#' 
#+ FigS5d
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
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/anticons_organism_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pd <- p$de_tool + 
  labs(y = "Proportion of anti-conservative\np histograms",
       x = "DE analysis tool")

#' 
#+ FigS5e
f <- anticons ~ de_tool + (1 | model)
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons_detool__1_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pe <- p$de_tool + 
  labs(y = "Proportion of anti-conservative p histograms",
       x = "DE analysis tool") + 
  theme(axis.title.y = element_blank())

#' 
#+ FigS5f
f <- anticons ~ de_tool + (de_tool | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons_detool__detool_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pf <- p$de_tool + 
  labs(y = "Proportion of anti-conservative p histograms",
       x = "DE analysis tool") + 
  theme(axis.title.y = element_blank())

#+ FigS5, fig.cap="Binomial logistic models for proportion of anti-conservative p histograms. A, simple model $anticons \\sim DEtool$. B, simple model $anticons \\sim DEtool$ fitted on complete data, N = 6267. C, model conditioned on year of GEO submission: $anticons \\sim year + DEtool$, N = 2109. D, model conditioned on studied organism (human/mouse/other): $anticons \\sim organism + DEtool$, N = 1733. E, varying intercept model $anticons \\sim DEtool + (1 | model)$ where 'model' stands for sequencing instrument model, N = 1718. F, varying intercept and slope model $anticons \\sim DEtool + (DEtool | model)$. B, C, E, and F y-axis labels are same as in A and D."
(pa + pb + pc) / (pd + pe + pf) + plot_annotation(tag_levels = "A")

#' 
#+ FigS6a
f <- pi0 ~ de_tool
family <- student()
data <- pvalues
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/pi0_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pa <- p$de_tool + 
  labs(x = "DE analysis tool")

#' 
#+ FigS6b
f <- pi0 ~ year + de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/pi0_year_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pb <- p$de_tool + 
  labs(x = "DE analysis tool")

#' 
#+ FigS6c
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == 9606 ~ "human",
    tax_id == 10090 ~ "mouse",
    TRUE ~ "other"
  ))
f <- pi0 ~ organism + de_tool
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/pi0_organism_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pc <- p$de_tool + 
  labs(x = "DE analysis tool")

#' 
#+ FigS6d
f <- pi0 ~ de_tool + (1 | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0_detool__1_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pd <- p$de_tool + 
  labs(x = "DE analysis tool")

#' 
#+ FigS6e
f <- pi0 ~ de_tool + (de_tool | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0_detool__detool_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pe <- p$de_tool + 
  labs(x = "DE analysis tool")

#+ FigS6, fig.cap="Robust (student-t likelihood) modeling of pi0. A, simple model $pi0 \\sim DEtool$ fitted on complete data, N = 1567. B, model conditioned on year of GEO submission: $pi0 \\sim year + DEtool$, N = 488. C, model conditioned on studied organism (human/mouse/other): $pi0 \\sim organism + DEtool$, N = 400. D, varying intercept model $pi0 \\sim DEtool + (1|model)$ where 'model' stands for sequencing instrument model, N = 396. E, varying intercept/slope model $pi0 \\sim DEtool + (DEtool | model)$, N = 396."
(pa + pb + pc) / (pd + pe + plot_spacer()) + plot_annotation(tag_levels = "A")

#' 
#+ FigS7, fig.cap="Modeling dependency of pi0 on sequencing instrument model: $pi0 \\sim (1 | model)$, N = 396."
f <- pi0 ~ (1 | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0__1_model.rds"))
mod %>%
  spread_draws(b_Intercept, r_model[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_model) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "pi0", y = "Sequencing instrument model") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS8, fig.cap="Modeling dependency of pi0 on library strategy: $pi0 \\sim (1| librarystrategy)$."
f <- pi0 ~ (1 | library_strategy)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0__1_librarystrategy.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_strategy[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_strategy) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "pi0", y = "Library strategy") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS9, fig.cap="Modeling dependency of pi0 on library selection: $pi0 \\sim (1| libraryselection)$, N = 396."
f <- pi0 ~ (1 | library_selection)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0__1_libraryselection.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_selection[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_selection) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "pi0", y = "Library selection") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+  FigS10, fig.cap="Modeling dependency of pi0 on library layout: $pi0 \\sim (1 | librarylayout)$, N = 396."
f <- pi0 ~ (1 | library_layout)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           iter = 4000, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0__1_librarylayout.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_layout[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_layout) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "pi0", y = "Library layout") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS11, fig.cap="Modeling dependency of pi0 on library layout: $pi0 \\sim (1 | library_source)$, N = 396."
f <- pi0 ~ (1 | library_source)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           iter = 4000, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/pi0__1_librarysource.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_source[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_source) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "pi0", y = "Library source") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS12, fig.cap="Modeling dependency of proportion of anti-conservative histograms on sequencing platform: $anticons \\sim (1  | model)$, N = 1718."
f <- anticons ~ (1 | model)
family <- bernoulli()
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons__1_model.rds"))
mod %>%
  spread_draws(b_Intercept, r_model[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_model) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative p histograms", y = "Sequencing instrument model") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS13, fig.cap="Modeling dependency of proportion of anti-conservative histograms on library strategy: $anticons \\sim (1  | librarystrategy)$, N = 1718."
f <- anticons ~ (1 | library_strategy)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons__1_librarystrategy.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_strategy[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_strategy) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative p histograms", y = "Library strategy") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS14, fig.cap="Modelling dependency of proportion of anti-conservative histograms on library selection: $anticons \\sim (1 | libraryselection)$, N = 1718."
f <- anticons ~ (1 | library_selection)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons__1_libraryselection.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_selection[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_selection) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative p histograms", y = "Library_selection") +
  scale_x_continuous(limits = c(0, 1))

#' 
#+ FigS15, fig.cap="Modeling dependency of proportion of anti-conservative histograms on library layout: $anticons \\sim (1 | librarylayout)$, N = 1718."
f <- anticons ~ (1 | library_layout)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("models/anticons__1_librarylayout.rds"))
mod %>%
  spread_draws(b_Intercept, r_library_layout[condition,]) %>%
  median_hdi(condition_mean = b_Intercept + r_library_layout) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative p histograms", y = "Library layout")
