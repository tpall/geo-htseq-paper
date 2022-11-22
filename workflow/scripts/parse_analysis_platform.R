library(tidyverse)
library(lubridate)
library(here)

document_summaries <- read_csv(here("results/data/document_summaries.csv"))
parsed_suppfiles_raw <- read_csv(here("results/data/parsed_suppfiles.csv"))
detools_from_pmc_raw <- read_csv(here("results/detools_from_pmc.csv"), col_types = "ccc")

acc_year <- document_summaries %>% 
  select(Accession, PDAT) %>% 
  mutate(year = year(PDAT))

parsed_suppfiles <- parsed_suppfiles_raw %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

get_var <- function(x) {
  NAs <- is.na(x)
  if (all(NAs)) {
    return("unknown")
  }
  return(names(x)[!NAs])
}

parsed_suppfiles_var <- parsed_suppfiles %>% 
  left_join(acc_year) %>% 
  select(Accession, id, Set, Type, Class, PDAT) %>% 
  filter(!is.na(Type)) %>% 
  pivot_wider(names_from = Type, values_from = Class) %>% 
  group_by(id, Set) %>% 
  mutate(var = get_var(c("logcpm"=logcpm, "rpkm"=rpkm, "fpkm"=fpkm, "basemean"=basemean, "aveexpr"=aveexpr))) %>%
  ungroup() %>% 
  mutate_at("Set", str_remove_all, '"') %>% 
  rowwise()

platform_regex <- "deseq2?|de(g|x)seq|rockhopper|cuff(diff|links)|edger|clc(bio)? ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth|partek|(nrsa|nascent rna seq)|median ratio norm|rmats|ballgown|biojupie|seurat|exdega"
normalise_name <- function(x) {
  x <- str_to_lower(x)
  x <- str_replace_all(x, "clc(bio)? ?genomics", "clc genomics")
  x <- str_replace_all(x, "(nrsa|nascent rna seq)", "nascent rna seq")
  x <- str_replace_all(x, "voom", "limma")
  return(x)
}

analysis_platform_from_varset <- parsed_suppfiles_var %>% 
  mutate(
    analysis_platform_from_set = case_when(
      str_detect(str_to_lower(Set), platform_regex) ~ str_extract(str_to_lower(Set), platform_regex), 
      str_detect(str_to_lower(Set), "tagwise") ~ "edger.tagwise",
      str_detect(str_to_lower(Set), "common") ~ "edger.common",
      str_detect(str_to_lower(Set), "baggerley") ~ "clc genomics",
      TRUE ~ NA_character_
      ),
    analysis_platform_from_var = case_when(
      var == "basemean" & str_detect(str_to_lower(Set), "pval(?!ue)") ~ "deseq",
      var == "basemean" & str_detect(str_to_lower(Set), "pvalue") ~ "deseq2",
      var == "logcpm" ~ "edger",
      var == "aveexpr" & str_detect(str_to_lower(Set), "p\\.value") & PDAT > "2014-01-01" ~ "limma",
      var == "fpkm" & str_detect(str_to_lower(Set), "p_value") ~ "cuffdiff",
      TRUE ~ NA_character_
      ),
    analysis_platform_from_filename = str_extract(str_to_lower(id), platform_regex)
    ) %>% 
  mutate_at(vars(starts_with("analysis_platform_from_")), normalise_name)

analysis_platform_from_summaries <- function(x) {
  x %>% 
    map(str_to_lower) %>% 
    map(str_extract_all, platform_regex) %>% 
    map(unlist) %>% 
    map(unique)
}

reorder <- function(x) {
  str_c(sort(unique(unlist(str_split(x, ";")))), collapse = ";")
}

document_summaries_analysis_platform <- document_summaries %>% 
  mutate(
    analysis_platform_from_sum = map(summary, analysis_platform_from_summaries)
  ) %>% 
  select(Accession, analysis_platform_from_sum) %>% 
  unnest(analysis_platform_from_sum) %>% 
  unnest(analysis_platform_from_sum) %>% 
  mutate_at(vars(starts_with("analysis_platform_from_")), normalise_name) %>% 
  group_by(Accession) %>% 
  summarise(analysis_platform_from_sum = str_c(unique(analysis_platform_from_sum), collapse = ";")) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(analysis_platform_from_sum = reorder(analysis_platform_from_sum)) %>% 
  ungroup()

accession_pubmedids <- document_summaries %>% 
  select(Accession, PubMedIds) %>% 
  drop_na() %>% 
  distinct() %>% 
  mutate(PubMedIds = str_split(PubMedIds, ";")) %>% 
  unnest(PubMedIds)

detools_from_pmc <- detools_from_pmc_raw %>% 
  left_join(accession_pubmedids) %>% 
  select(Accession, analysis_platform_from_pmc = detool)  %>% 
  mutate_at(vars(starts_with("analysis_platform_from_")), normalise_name) %>% 
  group_by(Accession) %>% 
  summarise(analysis_platform_from_pmc = str_c(unique(analysis_platform_from_pmc), collapse = ";")) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(analysis_platform_from_pmc = reorder(analysis_platform_from_pmc)) %>% 
  ungroup()

analysis_platform_from_varset %>%
  left_join(document_summaries_analysis_platform) %>% 
  left_join(detools_from_pmc) %>% 
  add_count(Accession, id, Set, sort = TRUE)

parse_platform <- function(...) {
  args <- list(...)
  if (all(is.na(args))) {
    return(NA_character_)
  }
  args <- args[!is.na(args)]
  args.untag <- str_remove_all(args, "\\.tagwise")
  tools <- unique(unlist(str_split(args.untag, ";")))
  if (all(tools %in% c("cuffdiff", "cufflinks"))) {
    tools <- "cuffdiff"
  }
  return(str_c(reorder(tools), collapse = ";"))
}

analysis_platform_summary <- analysis_platform_from_varset %>%
  left_join(document_summaries_analysis_platform) %>% 
  left_join(detools_from_pmc) %>% 
  rowwise() %>% 
  mutate(
    analysis_platform_parsed = parse_platform(analysis_platform_from_var, analysis_platform_from_set, analysis_platform_from_filename, analysis_platform_from_sum, analysis_platform_from_pmc),
    analysis_platform_parsed1 = case_when(
      str_detect(analysis_platform_parsed, "^;") ~ analysis_platform_parsed,
      analysis_platform_parsed == "deseq;deseq2" & analysis_platform_from_var == "deseq2" ~ "deseq2",
      analysis_platform_parsed == "deseq;deseq2" & analysis_platform_from_var == "deseq" ~ "deseq",
      analysis_platform_from_filename == "deseq" & analysis_platform_from_var == "deseq2" & is.na(analysis_platform_from_set) ~ "deseq2",
      analysis_platform_from_pmc == "ballgown;limma" ~ "ballgown",
      str_detect(analysis_platform_parsed, "clc genomics") & str_detect(analysis_platform_parsed, "edger") ~ "clc genomics-edger",
      str_detect(analysis_platform_parsed, str_remove(analysis_platform_from_set, "\\.(common|tagwise)")) & (is.na(analysis_platform_from_filename) | str_remove(analysis_platform_from_set, "\\.(common|tagwise)") == analysis_platform_from_filename) ~ str_remove(analysis_platform_from_set, "\\.(common|tagwise)"),
      str_detect(analysis_platform_parsed, analysis_platform_from_filename) & (is.na(analysis_platform_from_set) | analysis_platform_from_filename == str_remove(analysis_platform_from_set, "\\.(common|tagwise)")) ~ analysis_platform_from_filename,
      (analysis_platform_from_var == analysis_platform_from_pmc | analysis_platform_from_var == "cuffdiff" & analysis_platform_from_pmc == "cufflinks") ~ analysis_platform_from_var,
      str_detect(str_replace(analysis_platform_from_pmc, "cufflinks", "cuffdiff"), analysis_platform_from_var) ~ analysis_platform_from_var,
      (str_detect(analysis_platform_from_pmc, "edger") & analysis_platform_from_var == "limma" | analysis_platform_from_pmc == "edger;limma") ~ "edger-limma",
      str_detect(analysis_platform_from_pmc, "^cufflinks") & str_detect(analysis_platform_from_pmc, "^[^;]*;[^;]*$") ~ str_replace(analysis_platform_from_pmc, ";", "-"),
      !str_detect(analysis_platform_from_pmc, ";") ~ analysis_platform_from_pmc,
      analysis_platform_from_pmc == "deseq;deseq2" ~ "deseq2",
      TRUE ~ analysis_platform_parsed
    ),
    analysis_platform_final = case_when(
      str_detect(analysis_platform_parsed1, ";") | analysis_platform_parsed1 %in% c("ebseq", "noiseq") ~ NA_character_, 
      analysis_platform_from_pmc == "deseq" & analysis_platform_from_var == "cuffdiff" ~ "cuffdiff",
      analysis_platform_parsed1 == "cufflinks" ~ "cuffdiff",
      TRUE ~ analysis_platform_parsed1),
    ) 

analysis_platform_summary %>%
  select(Accession, id, Set, analysis_platform = analysis_platform_final) %>% 
  ungroup() %>% 
  write_csv(here("results/parse_analysis_platform_output.csv"))
