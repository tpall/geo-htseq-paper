if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}

library(readr)
library(dplyr)
library(stringr)
library(here)

spots_raw <- read_csv(here("results/data/spots.csv")) %>% 
  rename_all(str_to_lower)

probs <- problems(spots_raw)

if (nrow(probs) > 0) {
  
  spots_raw <- spots_raw %>% 
    mutate_all(as.character)
  
  #' Fix missing columns
  columns_off <- spots_raw[probs$row, ]
  colnames(columns_off) <- c(setdiff(colnames(columns_off), "sample_name"), "sample_name")
  
  columns_ok <- spots_raw[-probs$row, ]
  
  spots_col_fix <- bind_rows(columns_ok, columns_off)
  
} else {
  spots_col_fix <- spots_raw
}




#' Fix cases with missing platform.
#+
seq_platform <- spots_col_fix %>% 
  mutate_at(vars(platform, model), str_to_lower) %>% 
  mutate(model = case_when(
    is.na(model) & !is.na(platform) ~ platform,
    TRUE ~ model
  ),
  model = case_when(
    !str_detect(model, platform) ~ str_c(platform, model, sep = " "),
    TRUE ~ model
  ))

sequencing_metadata <- seq_platform %>% 
  mutate_at("spots", as.numeric) %>% 
  group_by(geo_accession, library_strategy, library_source, library_selection, library_layout, platform, model, tax_id) %>% 
  summarise_at("spots", list(reads = mean)) %>% 
  ungroup()

sequencing_metadata %>% 
  write_csv(here("results/sequencing_metadata.csv"))

sequencing_metadata_unique_platform <- sequencing_metadata %>% 
  group_by(Accession = geo_accession) %>% 
  add_count() %>% 
  ungroup() %>% 
  filter(n == 1) %>% 
  select(Accession, everything()) %>% 
  select(-geo_accession, -n)

write_csv(sequencing_metadata_unique_platform, here("results/sequencing_metadata_unique_platform.csv"))
