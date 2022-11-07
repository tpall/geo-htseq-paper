library(tidyverse)
library(here)

batches <- list.files("results", pattern = "^pmcid_", full.names = TRUE) %>% 
  map(read_csv, col_types = "c")

batches %>% 
  bind_rows() %>% 
  filter(!is.na(PMCID)) %>% 
  mutate_at("PMCID", str_remove, "^PMC") %>% 
  write_csv(here("results/pmcid.csv"))
