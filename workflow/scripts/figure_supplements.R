#' ---
#' title: Figure supplements
#' date: ""
#' author: ""
#' output:
#'    rmarkdown::pdf_document:
#'        number_sections: FALSE
#'        toc: FALSE
#' urlcolor: blue
#' header-includes:
#' - \usepackage[font=footnotesize,labelfont=bf]{caption}
#' - \DeclareCaptionLabelFormat{supplement}{S#2 Fig.}
#' - \captionsetup{labelformat=supplement,labelsep=quad}
#' bibliography: main/references.bib
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
library(here)
library(glue)
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
if (!dir.exists("results/models")) {
  message("Creating results/models dir..")
  dir.create("results/models", recursive = TRUE)
}

#+ data
pvalues <- read_csv(here("results/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_filtered <- read_csv(here("results/pvalues_filtered.csv")) %>% 
  rename(de_tool = analysis_platform)
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))
de_simulation_results <- read_csv(here("results/data/de_simulation_results.csv"))

#' 
#+ s1fig
number_of_pvalues <- pvalues_sample %>%
  select(Accession, id, hist, Class) %>% 
  rowwise() %>% 
  mutate(
    hist = sum(as.numeric(str_trim(str_extract_all(hist, "[0-9 ]+")[[1]]))),
    anticons = if_else(Class %in% c("anti-conservative", "uniform"), "anti-conservative", "all other classes")
  ) %>% 
  distinct() %>% 
  rename(n_pvalues = hist)

pa <- number_of_pvalues %>% 
  mutate(
    anticons = ifelse(anticons == "anti-conservative", "anti-conservative\nand uniform", anticons)) %>% 
  ggplot() +
  geom_histogram(aes(n_pvalues), binwidth = 0.1) +
  geom_vline(xintercept = median(number_of_pvalues$n_pvalues), linetype = "dashed") +
  facet_wrap(~anticons, scales = "free_y") +
  scale_x_log10() +
  labs(x = expression(Number~of~p~values~(log[10])), y = "Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

data <- number_of_pvalues %>% 
  mutate(log10_n_pvalues = log10(n_pvalues))
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
  select(anticons) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  mutate(anticons = ifelse(anticons == "anti-conservative", "anti-conservative\nand uniform", "all\nother\nclasses"))
pb <- draws %>% 
  ggplot(aes(anticons, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(Number~of~p~values~(log[10])), y = "Count") +
  theme(axis.title.x = element_blank())

p <- pa + pb + plot_layout(widths = c(2/3, 1/3)) + plot_annotation(tag_levels = "A")
ggsave(here("figures/S1Fig.tiff"), plot = p, width = 12, height = 8, units = "cm", dpi = 300)

number_of_pvalues_formatted <- prettyNum(nrow(number_of_pvalues), big.mark=',')
median_number_of_pvalues_formatted <- prettyNum(round(median(number_of_pvalues$n_pvalues)), big.mark=',')

fig_cap <- glue("__Reduced number of features in anti-conservative and uniform p value sets.__ 
                (__A__) P value set size distribution. Dashed line denotes the median ({median_number_of_pvalues_formatted}) number of features. From each GEO series only one random set was considered, N = {number_of_pvalues_formatted} p value sets. 
                (__B__) Robust linear modeling of number of features in anti-conservative and uniform vs. non-anti-conservative p value sets [log10_n_pvalues ~ anticons], student’s t likelihood, N = {number_of_pvalues_formatted}. 
                Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible region, respectively.
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/log_n_pvalues%20~%20anticons.rds."
                )

#+ fig.cap=fig_cap
p


#+ s2fig
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
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_year.rds"),
           file_refit = "on_change")
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions),
          line_args = list(color = "black"),
          plot = FALSE)

fig_cap <- glue("__The increasing proportion of anti-conservative histograms.__ 
                Binomial logistic model [{deparse(f)}], N = {prettyNum(summary(mod)$nobs, big.mark=',')}. 
                Lines denote best fit of linear model. Shaded area denotes 95% credible region.
                The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_year.rds.")

#+ fig.cap=fig_cap
p$year + 
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Year")
ggsave(here("figures/S2Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ s3fig
f <- anticons ~ year + (year | de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_year__year_detool.rds"),
           file_refit = "on_change"
)

draws <- data %>% 
  add_epred_draws(mod, ndraws = 100)

p <- draws %>% 
  ggplot(aes(year, .epred, group = 1)) +
  stat_lineribbon(point_interval = median_qi, .width = 0.95, fill = "gray65")

fig_cap <- glue("__All differential expression analysis tools are associated with temporally increasing anti-conservative p value histograms__, 
                two-level binomial logistic model [{deparse(f)}], 
                N = {prettyNum(summary(mod)$nobs, big.mark=',')}.
                Lines denote best fit of linear model. Shaded area denotes 95% credible region.
                The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_year__year_detool.rds.")

#+fig.cap=fig_cap
p + 
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Year") +
  facet_wrap(~ de_tool) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(here("figures/S3Fig.tiff"), height = 12, width = 18, dpi = 300, units = "cm")

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
  filter(!(model %in% c("454 gs flx", "454 gs flx titanium", "dnbseq-g400", "dnbseq-g50", "illumina miniseq", "ion s5", "helicos heliscope"))) %>% 
  add_epred_draws(mod, ndraws = 100)

p <- draws %>% 
  ggplot(aes(year, .epred, group = 1)) +
  stat_lineribbon(point_interval = median_qi, .width = 0.95, fill = "gray65") +
  facet_wrap(~ model, labeller = label_wrap_gen(width = 18))

fig_cap <- glue("__All sequencing instrument models are associated with temporally increasing anti-conservative p value histograms__,
                two-level binomial logistic model [{deparse(f)}], N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Only GEO submissions utilizing 
                single sequencing platform were used for model fitting. Lines denote best fit of linear model. Shaded area denotes 95% credible region. 
                The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_year__year_model.rds.")

#+ fig.height=8, fig.cap=fig_cap
p + 
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(y = "Proportion of anti-conservative p value histograms",
       x = "Year") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 8))
ggsave(here("figures/S4Fig.tiff"), height = 23, width = 18, dpi = 300, units = "cm")

#'
#+ s5fig
p <- pvalues_sample %>% 
  count(year, de_tool) %>% 
  add_count(year, name = "total", wt = n) %>% 
  mutate(Proportion = n / total) %>% 
  ggplot() +
  geom_line(aes(year, Proportion, color = de_tool), size = 1) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  labs(x = "Year") +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.8))

fig_cap <- glue("__No single differential expression analysis tool dominates the field.__ 
                Y-axis shows the proportion of analysis platforms, 
                x-axis shows publication year of GEO submission, N = {prettyNum(sum(p$data$n), big.mark=',')}.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S5Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

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
  select(de_tool) %>% 
  distinct() %>% 
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
  select(de_tool) %>% 
  distinct() %>% 
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
  select(de_tool, year) %>%
  filter(year == 2020) %>% 
  distinct() %>% 
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
  select(de_tool, organism) %>%
  filter(organism == "human") %>% 
  distinct() %>% 
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
  inner_join(sequencing_metadata)
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
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
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool_all.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_year_detool.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_organism_detool.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool__1_model.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool__detool_model.rds.')

#+ s6fig, fig.cap=fig_cap
(pa + pb + pc) / (pd + pe + pf) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))
ggsave(here("figures/S6Fig.tiff"), height = 14, width = 18, dpi = 300, units = "cm")

#' 
#+ s7figa
f <- pi0 ~ de_tool
family <- Beta()
data <- pvalues
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
  select(de_tool) %>% 
  distinct() %>% 
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
data <- pvalues_sample
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
  select(de_tool, year) %>%
  filter(year == 2020) %>% 
  distinct() %>% 
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
  ))
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
  select(de_tool, organism) %>%
  filter(organism == "human") %>% 
  distinct() %>% 
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
           data = data, 
           family = family, 
           prior = prior("student_t(3, 0, 1)", class = "b"),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0_detool__1_model.rds"),
           file_refit = "on_change")

draws <- as.data.frame(mod)
pd <- draws %>% 
  mutate_at(vars(starts_with("b_de_tool")), ~.x + draws$b_Intercept) %>% 
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

p <- (pa + pb + pc) / (pd + pe + plot_spacer()) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))
fig_cap <- glue("__DE analysis tool conditional effects from beta regression modeling of $\\pi_0$.__ 
                (__A__) Simple model [{paF}] fitted on complete data, N = {paN}. 
                (__B__) Model conditioned on year of GEO submission [{pbF}], N = {pbN}. 
                (__C__) Model conditioned on studied organism (human/mouse/other) [{pcF}], N = {pcN}. 
                (__D__) Varying intercept model [{pdF}] where 'model' stands for sequencing instrument model, N = {pdN}. 
                (__E__) Varying intercept/slope model [{peF}], N = {peN}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively.
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool_full_data.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_year_detool.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_organism_detool.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool__1_model.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool__detool_model.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S7Fig.tiff"), height = 14, width = 18, dpi = 300, units = "cm")

#' 
#+ s8fig
f <- pi0 ~ model
data <- data %>% 
  filter(!is.na(pi0), !is.na(model))
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
  select(model) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(model, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0])) +
  theme(axis.title.y = element_blank())

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on sequencing instrument model__ [{deparse(f)}], beta distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0__model.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S8Fig.tiff"), height = 12, width = 18, dpi = 300, units = "cm")

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
  select(library_strategy) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(library_strategy, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0]), y = "Library strategy")

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on library strategy__ [{deparse(f)}], beta distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0__librarystrategy.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S9Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

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
  select(library_selection) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(library_selection, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0]), y = "Library selection")

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on library selection__ [{deparse(f)}], beta distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0__libraryselection.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S10Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

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
  select(library_layout) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(library_layout, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = expression(pi[0]), y = "Library layout")

fig_cap <- glue("__Modeling dependency of $\\pi_0$ on library layout__ [{deparse(f)}], beta distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0__librarylayout.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S11Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ s12fig
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  filter(!is.na(model))
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
  select(model) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(model, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms") +
  theme(axis.title.y = element_blank())

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on sequencing platform__ [{deparse(f)}], bernoulli distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons__model.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S12Fig.tiff"), height = 12, width = 18, dpi = 300, units = "cm")

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
  select(library_strategy) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(library_strategy, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms", y = "Library strategy")

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on library strategy__ [{deparse(f)}], bernoulli distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons__librarystrategy.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S13Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

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
  select(library_selection) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(library_selection, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms", y = "Library selection")

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on library selection__ [{deparse(f)}], bernoulli distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons__libraryselection.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S14Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

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
  select(library_layout) %>% 
  distinct() %>% 
  add_epred_draws(mod) %>% 
  ggplot(aes(.epred, reorder(library_layout, .epred, median))) +
  stat_pointinterval(point_size = 1) +
  labs(x = "Proportion of anti-conservative histograms", y = "Library layout")

fig_cap <- glue("__Modeling dependency of proportion of anti-conservative histograms on library layout__ [{deparse(f)}], bernoulli distribution, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Thick and thin lines denote 66% and 95% credible interval, respectively. The model object related to figure can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons__librarylayout.rds.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/S15Fig.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

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
  select(de_tool) %>% 
  distinct() %>% 
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
  select(de_tool) %>% 
  distinct() %>% 
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
  select(de_tool, year) %>% 
  filter(year == 2020) %>% 
  distinct() %>% 
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
  select(de_tool, organism) %>% 
  filter(organism == "human") %>% 
  distinct() %>% 
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  # scale_y_continuous(limits = c(0, 1)) +
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
  # scale_y_continuous(limits = c(0, 1)) +
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
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool_filtered.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool_all_filtered.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_year_detool_filtered.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_organism_detool_filtered.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool__1_model_filtered.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/anticons_detool__detool_model_filtered.rds.')

#+ s16fig
p16 <- (pa + pb + pc) / (pd + pe + pf) +  
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))
ggsave(here("figures/S16Fig.tiff"), plot = p16, height = 15, width = 18, dpi = 300, units = "cm")

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

draws_17aa <- pvalues_filtered_sample %>% 
  select(de_tool) %>% 
  distinct() %>% 
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
es_pi0_prop <- pvalues_sample %>% 
  select(de_tool) %>% 
  distinct() %>% 
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

draws_17a <- data %>% 
  select(de_tool) %>% 
  distinct() %>% 
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

pb <- data %>% 
  select(de_tool, year) %>% 
  filter(year == 2020) %>% 
  distinct() %>% 
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

pc <- data %>% 
  select(de_tool, organism) %>% 
  filter(organism == "human") %>% 
  distinct() %>% 
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
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
  select(cuffdiff = b_Intercept, starts_with("b_de_tool")) %>% 
  rename_all(str_remove, "b_de_tool") %>% 
  pivot_longer(cols = colnames(.), names_to = "de_tool") %>% 
  mutate_at("value", inv_logit_scaled) %>% 
  ggplot(aes(de_tool, value)) +
  stat_pointinterval(point_size = 1) +
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
                The model object related to panel A can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool_sample_filtered.rds. 
                The model object related to panel B can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool_full_data_filtered.rds. 
                The model object related to panel C can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_year_detool_filtered.rds. 
                The model object related to panel D can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_organism_detool_filtered.rds. 
                The model object related to panel E can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool__1_model_filtered.rds. 
                The model object related to panel F can be downloaded from https://gin.g-node.org/tpall/geo-htseq-paper/raw/26619a4b74aa3781ac6a244edcc24e0ad6eb064b/models/pi0_detool__detool_model_filtered.rds.")

#+ 
ggsave(here("figures/S17Fig.tiff"), plot = p17, height = 15, width = 18, dpi = 300, units = "cm")


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
ggsave(here("figures/S18Fig.tiff"), plot = p, width = 12, height = 8, units = "cm", dpi = 300)

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


