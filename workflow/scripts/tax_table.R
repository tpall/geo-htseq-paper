library(tidyverse)
library(taxize)
library(here)

sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))
taxize_class <- classification(na.omit(unique(sequencing_metadata$tax_id)), db = "ncbi")

tax_table <- taxize_class %>% 
  map(filter, rank %in% c("superkingdom", "kingdom", "phylum", "order", "family", "genus", "species")) %>% 
  map(select, c("name", "rank")) %>% 
  map(pivot_wider, names_from = "rank", values_from = "name") %>% 
  bind_rows(.id = "tax_id")
tax_table %>% 
  write_csv(here("results/tax_table.csv"))
