#' ---
#' title: Supplementary figures
#' date: ""
#' author: ""
#' output:
#'    bookdown::pdf_document2:
#'        number_sections: FALSE
#'        toc: FALSE
#' ---

#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, message = FALSE, comment = FALSE, warning = FALSE)


#+ libs
library(readr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(cowplot)
library(patchwork)
library(brms)
library(here)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

#'
#+ params
chains <- 4
cores <- chains
refresh = 0

#+ data
pvalues <- read_csv(here("output/pvalues.csv"))
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
                              conditions = conditions, 
                              re_formula = NULL),
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
de_tool <- unique(data$de_tool)
conditions <- data.frame(de_tool = de_tool, row.names = de_tool)
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
model <- na.omit(unique(data$model))
conditions <- data.frame(model, row.names = model)
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
                              effects = "de_tool",
                              re_formula = NULL),
          plot = FALSE)
pa <- p$de_tool + 
  labs(y = "Proportion of anti-conservative\np histograms",
       x = "DE analysis tool")

#'
#+ FigS5b
data <- pvalues %>% 
  rename(de_tool = analysis_platform)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("models/anticons_detool_all.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool",
                              re_formula = NULL),
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
                              effects = "de_tool",
                              re_formula = NULL),
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
                              effects = "de_tool",
                              re_formula = NULL),
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
                              effects = "de_tool",
                              re_formula = NULL),
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
                              effects = "de_tool",
                              re_formula = NULL),
          plot = FALSE)
pf <- p$de_tool + 
  labs(y = "Proportion of anti-conservative p histograms",
       x = "DE analysis tool") + 
  theme(axis.title.y = element_blank())

#+ FigS5, fig.cap="Binomial logistic models for proportion of anti-conservative p histograms. A, simple model $anticons \\sim DEtool$. B, simple model $anticons \\sim DEtool$ fitted on complete data, N = 6267. C, model conditioned on year of GEO submission: $anticons \\sim year + DEtool$, N = 2109. D, model conditioned on studied organism (human/mouse/other): $anticons \\sim organism + DEtool$, N = 1733. E, varying intercept model $anticons \\sim DEtool + (1 | model)$ where 'model' stands for sequencing instrument model, N = 1718. F, varying intercept and slope model $anticons \\sim DEtool + (DEtool | model)$. B, C, E, and F y-axis labels are same as in A and D."
(pa + pb + pc) / (pd + pe + pf) + plot_annotation(tag_levels = "A")


#' ## Figure S6. Robust (student-t likelihood) modeling of pi0. 
#' A. Simple model pi0~platform fitted on complete data. 
#' 
#+
f <- pi0 ~ platform
family <- student()
data <- "all_pvalues"

#' C. Model conditioned on year of GEO submission: pi0~year + platform. 
#' 
#+
f <- pi0 ~ year + platform
data <- "pvalues_sample"

#' D. Model conditioned on studied organism (human/mouse/other): pi0~organism + platform. 
#' 
#+
f <- pi0 ~ organism + platform

#' E. Varying intercept model pi0~platform + (1|model) where “model” stands for sequencing platform. 
#' 
#+
f <- pi0 ~ platform + (1 | model)

#' F. Varying intercept/slope model pi0~platform + (platform|model).
#' 
#+
f <- pi0 ~ platform + (platform | model)

#' ## Figure S7. Modeling dependency of pi0 on sequencing platform: pi0~(1|model).
#' 
#+
f <- pi0 ~ (1 | model)

#' ## Figure S8. Modeling dependency of pi0 on library strategy: pi0~(1| library_strategy).
#' 
#+
f <- pi0 ~ (1 | library_strategy)

#' Figure S9. Modeling dependency of pi0 on library selection: pi0~(1| library_selection).
#' 
#+
f <- pi0 ~ (1 | library_selection)

#' ## Figure S10. Modeling dependency of pi0 on library layout: pi0~(1| library_layout).
#' 
#+
f <- pi0 ~ (1 | library_layout)

#' ## Figure S11. Modeling dependency of pi0 on library layout: pi0~(1| library_source).
#' 
#+
f <- pi0 ~ (1 | library_source)

#' ## Figure S12. Modeling dependency of proportion of anti-conservative histograms on 
#' sequencing platform: anticons~(1|model).
#' 
#+
f <- anticons ~ (1 | model)
family <- bernoulli()


#' ## Figure S13. Modeling dependency of proportion of anti-conservative histograms on 
#' library strategy: anticons~(1| library_strategy).
#' 
#+
f <- anticons ~ (1 | library_strategy)


#' ## Figure S14. Modelling dependency of proportion of anti-conservative histograms on 
#' library selection: anticons~(1| library_selection).
#' 
#+
f <- anticons ~ (1 | library_selection)

#' ## Figure S15. Modeling dependency of proportion of anti-conservative histograms on 
#' library layout: anticons~(1| library_layout).
#' 
#+
f <- anticons ~ (1 | library_layout)


