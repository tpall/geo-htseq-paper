if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}

library(stats) # load first otherwise masks dplyr filter
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(lubridate)
library(here)

#' #' Number of unique GEOs
#' #+
document_summaries <- read_csv(here("results/data/document_summaries.csv"))

if (!exists("snakemake")) {
  document_summaries %>%
    filter(PDAT<="2020-12-31") %>%
    pull(Accession) %>%
    n_distinct()
}

acc_year <- document_summaries %>%
  select(Accession, PDAT) %>%
  mutate(year = year(PDAT))

#' Number of sets with p-values
#+
col_types <- cols(
  id = col_character(),
  Type = col_character(),
  Class = col_character(),
  Conversion = col_character(),
  pi0 = col_double(),
  FDR_pval = col_double(),
  hist = col_character(),
  note = col_character(),
  Set = col_character()
)
parsed_suppfiles_raw <- read_csv(here("results/data/parsed_suppfiles.csv"), col_types = col_types)
parsed_suppfiles <- parsed_suppfiles_raw %>% 
  filter(!str_detect(id, "_RAW.tar")) %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

#' Parse analysis platform
get_var <- function(x) {
  NAs <- is.na(x)
  if (all(NAs)) {
    return("unknown")
  }
  return(names(x)[!NAs])
}

parsed_suppfiles_filtered <- parsed_suppfiles %>% 
  filter(Type != "raw")

pi0 <- parsed_suppfiles_filtered %>% 
  select(Accession, id, Set, pi0) %>% 
  na.omit()

hist <- parsed_suppfiles_filtered %>% 
  select(Accession, id, Set, FDR_pval, hist)

pvalues <- parsed_suppfiles_filtered %>% 
  select(Accession, id, Set, Type, Class) %>%
  mutate( # parsing analysis platform using expression level variable name
    analysis_platform_from_expression = case_when(
      Type == "basemean" ~ "deseq",
      Type == "aveexpr" ~ "limma",
      Type == "logcpm" ~ "edger",
      Type == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
      TRUE ~ "unknown"
    ), # parsing analysis platform using file names
    analysis_platform_from_filename = case_when(
      str_detect(str_to_lower(id), "deseq") ~ "deseq",
      str_detect(str_to_lower(id), "edger") ~ "edger",
      str_detect(str_to_lower(id), "limma") ~ "limma",
      str_detect(str_to_lower(id), "cuff") ~ "cuffdiff",
      TRUE ~ "unknown"
    ), # assigning analysis platform using file names and expression level variable name
    analysis_platform = case_when(
      analysis_platform_from_expression == analysis_platform_from_filename ~ analysis_platform_from_expression,
      analysis_platform_from_filename == "unknown" ~ analysis_platform_from_expression,
      TRUE ~ analysis_platform_from_filename
    )
  ) %>%
  select(Accession, id, Set, Class, analysis_platform) %>% 
  left_join(acc_year) %>% 
  left_join(pi0) %>% 
  left_join(hist) %>% 
  select(-PDAT) %>% 
  mutate(anticons = case_when(
    Class == "anti-conservative" ~ 1,
    TRUE ~ 0
  ))
write_csv(pvalues, here("results/pvalues_filtered.csv"))

#' P values
#+

#' Number of distinct GEO-s with p values
if (!exists("snakemake")) {
  n_distinct(pvalues$Accession)
}

#' Sample 1 p value set per accession
set.seed(11)
pvalues_sample <- pvalues %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  ungroup()
write_csv(pvalues_sample, here("results/pvalues_filtered_sample.csv"))

#' Save table ids for figure
#+
pvalues_sample %>% 
  select(id, Set) %>% 
  write_csv(here("results/pvalues_filtered_acc.csv"))

if (!exists("snakemake")) {
  pvalues_sample %>% 
    count(Class) %>% 
    summarise(Class, 
              n,
              p = signif(n / sum(n), 2))
}
