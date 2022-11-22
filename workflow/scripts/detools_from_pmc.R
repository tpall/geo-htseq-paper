library(tidyverse)
library(here)

document_summaries <- read_csv(here("results/data/document_summaries.csv")) %>% 
  select(Accession, PubMedIds) %>% 
  distinct() %>% 
  drop_na()

pmcid <- read_csv(here("results/pmcid.csv"), col_types = "cc") %>% 
  rename(PubMedIds = PMID) %>% 
  select(PubMedIds, PMCID) %>% 
  distinct()

detools_from_pmc_raw <- list.files("results", pattern = "detools_pmcid", full.names = TRUE) %>% 
  map(read_csv, col_types = "cccc") %>% 
  bind_rows() %>% 
  left_join(pmcid) %>% 
  filter(!is.na(detool))

detools_from_pmc_raw2 <- list.files("results", pattern = "detools_from_pubmed2", full.names = TRUE) %>% 
  map(read_csv, col_types = "cccc") %>% 
  bind_rows() %>% 
  left_join(pmcid) %>% 
  filter(!is.na(detool))

detools_from_pmc_raw3 <- list.files("results", pattern = "detools_from_pubmed3", full.names = TRUE) %>% 
  map(read_csv, col_types = "cccc") %>% 
  bind_rows() %>% 
  left_join(pmcid) %>% 
  filter(!is.na(detool))



detools_from_series_matrix <- read_csv(here("results/detools_from_series_matrix.csv"), col_types = "ccc") %>% 
  drop_na()

detools_pmcid_4 <- list.files("results", pattern = "detools_pmcid_4", full.names = TRUE) %>% 
  map(read_csv, col_types = "cccc") %>% 
  bind_rows() %>% 
  left_join(pmcid) %>% 
  filter(!is.na(detool))

detools_from_series_matrix_and_pmc_updated_regex <- read_csv(here("results/detools_from_series_matrix_and_pmc_updated_regex.csv"), col_types = "cc") %>% 
  left_join(pmcid) %>% 
  filter(!is.na(detool))



detools_from_pmc <- detools_from_pmc_raw %>% 
  bind_rows(detools_from_pmc_raw2) %>% 
  bind_rows(detools_from_pmc_raw3) %>% 
  bind_rows(detools_from_series_matrix) %>% 
  bind_rows(detools_pmcid_4) %>%
  bind_rows(detools_from_series_matrix_and_pmc_updated_regex) %>% 
  select(PubMedIds, PMCID, detool) %>% 
  group_by(PubMedIds, PMCID) %>% 
  summarise(detool = str_c(unique(str_split(str_c(detool, collapse = ";"), ";")[[1]]), collapse = ";")) %>% 
  ungroup() %>% 
  distinct()

detools_from_pmc %>% 
  write_csv(here("results/detools_from_pmc.csv"))

# detools_from_pmc %>% 
#   count(detool, sort = TRUE)

# d <- list.files("results", pattern = "detools_pmcid", full.names = TRUE) %>% 
#   map(read_csv, col_types = "cccc") %>% 
#   bind_rows() %>% 
#   left_join(pmcid)
# 
# 
# d %>% 
#   filter(str_detect(context, "pmc") | is.na(context)) %>% 
#   anti_join(detools_from_pmc_raw2 %>% select(PubMedIds, PMCID)) %>% 
#   filter(!is.na(context)) %>% 
#   write_csv(here("results/pubmedids_full_restricted.csv"))
