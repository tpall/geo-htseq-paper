library("tidyverse")
library("sparkline")
library("formattable")
library("here")

#' Parse GEO to sample accession mapping. 
#+
document_summaries <- read_csv("data/document_summaries.csv")
sample_acc_map <- document_summaries %>% 
  mutate(sample_acc = str_extract_all(Samples, "GSM\\d+")) %>% 
  select(accession = Accession, sample_acc)
sample_acc <- sample_acc_map$sample_acc
accession <- rep(sample_acc_map$accession, map_int(sample_acc, length))
acc_map_unnested <- tibble(accession, sample_acc = flatten_chr(sample_acc))

#' Import read runs data.
#+
pvals_wide <- read_csv("output/pvals_wide.csv")
spots_imported <- read_csv("data/spots.csv", col_types = cols(tax_id = col_character()))

#' Fix as much names as possible, parse title, strategy, and organism.
#' 'Name' field in older series is not structured by ";" separator.
#' We have ~417 series where we can't parse these fields from name.
#+ 
spots <- spots_imported %>% 
  rename_all(str_to_lower) %>% 
  filter(is.na(err)) %>% 
  select(-err)

#' Summarize missing raw data metadata error messages.
#+
spots_err <- spots_imported %>% 
  filter(!is.na(err)) %>% 
  select(geo_accession, err) %>% 
  mutate_at("err", str_replace, "SRA Experiment SRX\\d+ is not public", "SRA Experiment is not public")

spots_err %>% 
  count(err)

platforms <- spots %>% 
  select(platform, model) %>% 
  distinct() %>% 
  filter(!is.na(model))

spots <- spots %>% 
  mutate(model = case_when(
    is.na(model) ~ platform,
    TRUE ~ model
  ),
  platform = case_when(
    model == platform ~ NA_character_,
    TRUE ~ platform
  )) %>% 
  full_join(platforms, .) %>% 
  select(colnames(spots))

spots %>% 
  pull(library_layout) %>% 
  unique()

offset <- spots %>% 
  filter(!(library_layout %in% c("PAIRED", "SINGLE")))

spots %>% 
  filter((library_layout %in% c("PAIRED", "SINGLE"))) %>% 
  select(starts_with("library"), platform, model) %>% 
  map(unique)

#' Summarize merged spots data
#+
spot_summary <- spots %>% 
  group_by(accession = geo_accession, tax_id, organism, library_strategy, library_source, library_selection, library_layout, platform, model) %>% 
  summarise_at(vars(spots, bases), median)

#' Merge p values with spots data.
#+
pvals_spots <- pvals_wide %>% 
  select(accession, id, PDAT, starts_with("Class"), starts_with("pi0")) %>% 
  left_join(spot_summary)

pvals_spots %>% 
  group_by(accession) %>% 
  filter(all(is.na(spots))) %>% 
  ggplot() +
  geom_bar(aes(PDAT))

