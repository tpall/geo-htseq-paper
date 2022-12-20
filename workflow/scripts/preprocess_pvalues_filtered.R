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

parsed_suppfiles_filtered <- parsed_suppfiles %>% 
  filter(Type != "raw", !is.na(Type))

analysis_platf <- read_csv(here("results/pvalues.csv")) %>% 
  select(Accession, id, Set, analysis_platform)

acc_year <- read_csv(here("results/data/document_summaries.csv")) %>% 
  select(Accession, PDAT) %>% 
  mutate(year = year(PDAT))

pvalues <- parsed_suppfiles_filtered %>% 
  left_join(analysis_platf) %>% 
  left_join(acc_year) %>% 
  select(-PDAT, -note) %>% 
  mutate(anticons = case_when(
    Class == "anti-conservative" ~ 1,
    TRUE ~ 0
  ))
write_csv(pvalues, here("results/pvalues_filtered.csv"))
