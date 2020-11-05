
#+ libs
library(readr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(cowplot)
library(patchwork)
library(viridisLite)
library(brms)
library(tidybayes)
library(here)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

#+ params
chains <- 4
refresh <- 1
rstan_options(javascript = FALSE)

#+ data
conformity_acc <- read_csv(here("output/conformity_acc.csv"))
pvalues <- read_csv(here("output/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample <- read_csv(here("output/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
sequencing_metadata <- read_csv(here("output/sequencing_metadata_unique_platform.csv"))

#'
#+ fig1
f <- conforms ~ year
family <- bernoulli()
data <- conformity_acc
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           refresh = refresh,
           file = here("models/conforms_year.rds"))
p <- plot(conditional_effects(mod), plot = FALSE)$year
p + 
  geom_smooth(color = "black") +
  labs(x = "Year", y = "Proportion of submissions conforming\nwith GEO submission guidelines") +
  scale_x_continuous(breaks = seq(2006, 2019, by = 2)) +
  scale_y_continuous(limits = c(0, 1))
ggsave(here("plots/conforming_per_year.png"), height = 7, width = 11, dpi = 300, units = "cm")

#'
#+ fig3
f <- Class ~ year + (year | de_tool)
family <- categorical()
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = chains, 
           refresh = refresh,
           iter = 3000,
           file = here("models/Class_year__year_detool.rds"))
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
  labs(y = "Proportion",
       x = "Year")

#'
#'
#+
f <- Class ~ de_tool
family <- categorical()
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = chains, 
           refresh = refresh,
           iter = 3000,
           file = here("models/Class_detool.rds"))
p <- plot(conditional_effects(mod, 
                              categorical = TRUE, 
                              effects = "de_tool"), 
          plot = FALSE)
p3b <- p$`de_tool:cats__` + theme(axis.title.x = element_blank())

#'
#'
#+
f <- pi0 ~ de_tool
family <- student()
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = chains, 
           refresh = refresh,
           file = here("models/pi0_detool.rds"))
p <- plot(conditional_effects(mod, 
                              conditions = conditions,
                              re_formula = NA),
          plot = FALSE)
p4b <- p$de_tool +
  theme(axis.title.x = element_blank())

#+
f <- pi0 ~ year + (year | de_tool)
prior <- c(prior(normal(0, 0.5), class = b), prior(lkj(3), class = cor), prior(student_t(5, 0, 0.5), class = sigma))
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = 1, 
           cores = 1, 
           refresh = refresh,
           prior = prior,
           iter = 3000,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
           file = here("models/pi0_year__year_detool.rds"))

conditions <- make_conditions(data, vars = "de_tool")
rownames(conditions) <- conditions$de_tool
p <- plot(conditional_effects(mod, 
                              effects = "year",
                              conditions = conditions,
                              re_formula = NULL),
          plot = FALSE)
p4c <- p$year + 
  facet_wrap(~ de_tool) +
  scale_x_continuous(breaks = seq(2009, 2019, 2))


p4a <- pvalues_sample %>% 
  filter(Class == "anti-conservative") %>% 
  ggplot() + 
  geom_histogram(aes(pi0), color = "white", binwidth = 0.1)

#+ Fig4, fig.cap=""
p4 <- (p4a + p4b) / p4c + 
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "A")
ggsave(here("plots/figure_4.png"), p4)
