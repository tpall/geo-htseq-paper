library("tidyverse")
library("here")

#' ## P value stats
#' Import dataset
imported_raw <- read_csv(here("data/parsed_suppfiles.csv"))
imported <- imported_raw %>% 
  mutate(accession = str_extract(id, "GSE\\d+")) %>% 
  select(accession, everything())

#' ## Experiment metadata
#' Importing experiment metadata (dates and etc) 
document_summaries <- read_csv(here("data/document_summaries.csv"))
range(document_summaries$PDAT)


#' Sequencing metadata
#+
spots_raw <- read_csv(here("data/spots.csv"), col_types = cols(tax_id = col_character())) %>% 
  rename_all(str_to_lower)
spots <- spots_raw %>% 
  group_by(geo_accession, library_strategy, library_source, library_selection, library_layout, platform, model, organism) %>% 
  summarise_at("spots", list(reads = mean)) %>% 
  ungroup()

#' Fix cases with missing platform.
#+
seq_platform <- spots %>% 
  mutate_at(vars(platform, model), str_to_lower) %>% 
  mutate(model = case_when(
    is.na(model) & !is.na(platform) ~ platform,
    TRUE ~ model
  ),
  model = case_when(
    !str_detect(model, platform) ~ str_c(platform, model, sep = " "),
    TRUE ~ model
  ))

seq_platform %>% 
  filter(library_strategy == organism)

seq_platform %>% 
  filter(organism == "RNA-Seq")

seq_platform %>% 
  select(geo_accession, model, organism) %>% 
  distinct()

#' Number of unique GEO sets during period.
#+
n_distinct(document_summaries$Accession)

#' Number of GEO sets with supplementary files that were imported.
#+
n_distinct(imported$accession)

#' Importing publications metadata
publications <- read_csv(
  here("data/publications.csv"),
  col_types = str_c(rep("c", 26), collapse = "")
  ) %>% 
  rename(PubMedIds = Id)

#' Single-cell experiments
single_cell <- read_csv(here("data/single-cell.csv"))

#' Citations
citations <- read_csv(
  here("data/scopus_citedbycount.csv"), 
  col_types = cols(PubMedIds = col_character())
  )

#' Merging publications with citations
pubs <- publications %>% 
  left_join(citations) %>% 
  select(PubMedIds, PubDate, Source, FullJournalName, ISSN, ESSN, citations) %>% 
  mutate(ISSN = case_when(
    is.na(ISSN) ~ ESSN,
    TRUE ~ ISSN
  ))

#' ## Updating P values with some metadata
pvals_wide <- document_summaries %>% 
  select(accession = Accession, PDAT, taxon, PlatformTitle) %>% 
  inner_join(wide) %>% 
  mutate(single_cell = accession %in% single_cell$Accession) %>% 
  mutate(platform = case_when(
    filter_var == "basemean" ~ "deseq",
    filter_var == "aveexpr" ~ "limma",
    filter_var == "logcpm" ~ "edger",
    filter_var == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
    TRUE ~ "unknown"
  ))
write_csv(pvals_wide, "output/pvals_wide.csv")
range(pvals_wide$PDAT)

pvals_wide %>% 
  count(platform)

pvals_wide %>% 
  count(taxon) %>% 
  arrange(desc(n))
