#' ---
#' title: Figure supplements of filtered p values
#' date: ""
#' author: ""
#' output:
#'    rmarkdown::pdf_document:
#'        number_sections: FALSE
#'        toc: FALSE
#' header-includes:
#' - \usepackage{caption}
#' - \captionsetup[figure]{labelformat=empty}
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
old <- theme_set(theme_cowplot(font_size = 12, font_family = "Helvetica"))

#'
#+ params
is_ci <- function() {
  "CI" %in% Sys.getenv()
}
chains <- ifelse(is_ci(), 1, 4)
cores <- chains
refresh = 0
rstan_options(auto_write = TRUE, javascript = FALSE)
if (!dir.exists("results/models")) {
    message("Creating results/models dir..")
    dir.create("results/models", recursive = TRUE)
}

#+ data
col_types <- cols(
  Accession = col_character(),
  id = col_character(),
  Set = col_character(),
  Class = col_character(),
  analysis_platform = col_character(),
  year = col_double(),
  pi0 = col_double(),
  FDR_pval = col_double(),
  hist = col_character(),
  anticons = col_double()
)
pvalues <- read_csv(here("results/pvalues_filtered.csv"), col_types = col_types) %>% 
  rename(de_tool = analysis_platform)
col_types <- cols(
  Accession = col_character(),
  id = col_character(),
  Set = col_character(),
  Class = col_character(),
  analysis_platform = col_character(),
  year = col_double(),
  pi0 = col_double(),
  FDR_pval = col_double(),
  hist = col_character(),
  anticons = col_double()
)
pvalues_sample <- read_csv(here("results/pvalues_filtered_sample.csv"), col_types = col_types) %>% 
  rename(de_tool = analysis_platform)

col_types <- cols(
  Accession = col_character(),
  library_strategy = col_character(),
  library_source = col_character(),
  library_selection = col_character(),
  library_layout = col_character(),
  platform = col_character(),
  model = col_character(),
  tax_id = col_double(),
  reads = col_double()
)
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"), col_types = col_types)


#+ fig2supp1
qc_threshold <- function(n_pvals, bins, fdr) {
  qbinom(1 - 1 / bins * fdr, n_pvals, 1 / bins)
}
qc_thr <- qc_threshold(n_pvals = 20000, bins = 40, fdr = 0.05)

#'
#+ Fig5-Supp1A
y_title <- "Prop. anti-cons."
data <- pvalues
f <- anticons ~ de_tool
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_detool_all_filtered.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pa <- p$de_tool + 
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

#'
#+ Fig5-Supp1b
f <- pi0 ~ de_tool
family <- student()
data <- pvalues
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0_detool_full_data_filtered.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pb <- p$de_tool + 
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)


#' 
#+ Fig5-Supp1
for (p in list(pa, pb)) {
  p$layers[[1]]$aes_params$size <- 1
}
p <- pa + pb +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

fig_cap <- glue("Figure 5--figure supplement 1. Removal of low-count genes from p value sets. 
                A, Binomial logistic model for proportion of anti-conservative p value histograms, *{pbF}* model was fitted on complete data, N = {pbN}.
                B, Robust (student's t likelihood) modeling of $\\pi_0$, *{paF}* model was fitted on complete data, N = {paN}. 
                Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_5_figure_supplement_1.tiff"), height = 7, width = 12, dpi = 300, units = "cm")
