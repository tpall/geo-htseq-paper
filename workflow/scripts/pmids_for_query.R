library(tidyverse)
library(here)

document_summaries <- read_csv(here("results/data/document_summaries.csv"))
parsed_suppfiles_raw <- read_csv(here("results/data/parsed_suppfiles.csv"))
parsed_suppfiles <- parsed_suppfiles_raw %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything()) %>% 
  filter(!is.na(Set))

pmids_for_query <- document_summaries %>% 
  select(Accession, PubMedIds) %>% 
  drop_na() %>% 
  inner_join(parsed_suppfiles %>% select(Accession) %>% distinct()) %>% 
  mutate(PubMedIds = str_split(PubMedIds, ";")) %>% 
  unnest(PubMedIds) %>% 
  select(PubMedIds) %>% 
  distinct()

pmids_for_query %>% 
  write_csv(here("results/pmids_for_query.csv"))
