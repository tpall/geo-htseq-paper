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
  
  spots <- bind_rows(columns_ok, columns_off)
} else {
  spots <- spots_raw
}

#' Fix cases with missing platform.
#+
seq_platform <- spots %>% 
  mutate_at(vars(platform, model), str_to_lower) %>% 
  mutate(
    model = ifelse(is.na(model) & platform == "bgiseq-500", platform, model),
    platform = ifelse(platform == "bgiseq-500" & model == "bgiseq-500", "bgiseq", platform),
    model = ifelse(is.na(model) & platform == "hiseq x ten", platform, model),
    platform = ifelse(platform == "hiseq x ten" & model == "hiseq x ten", "illumina", platform),
    model = ifelse(is.na(model) & str_detect(platform, "illumina"), platform, model),
    platform = ifelse(str_detect(model, "illumina"), "illumina", platform),
    )

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
