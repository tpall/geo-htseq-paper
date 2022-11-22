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

#' Number of unique GEOs
#+
document_summaries <- read_csv(here("results/data/document_summaries.csv"))
if (!exists("snakemake")) {
  document_summaries %>% 
    filter(PDAT <= "2020-12-31") %>% 
    pull(Accession) %>% 
    n_distinct()
}

if (!exists("snakemake")) {
  document_summaries %>% 
    filter(PDAT <= "2020-12-31") %>% 
    group_by(year(PDAT)) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(cumsum = cumsum(n),
           perc = n / cumsum)
}

#' Supplementary files. We discard all 'filelist.txt' files.
#+
suppfilenames <- read_lines(here("results/data/suppfilenames.txt")) %>% 
  tibble(suppfilenames = .) %>% 
  filter(
    str_detect(suppfilenames, "suppl"),
    !str_detect(suppfilenames, "suppl/filelist.txt")
  )

acc_year <- document_summaries %>% 
  select(Accession, PDAT) %>% 
  mutate(year = year(PDAT))

#' Keep only subset from our time frame. Recheck inner_join!
conformity <- suppfilenames %>% 
  mutate(Accession = str_to_upper(str_extract(suppfilenames, "GS[Ee]\\d+"))) %>% 
  full_join(acc_year) %>%
  mutate(
    conforms = !(is.na(suppfilenames) | str_detect(str_to_lower(suppfilenames), c("readme", "\\.[bs]am", "\\.bed") %>% str_c("(\\.gz)?", collapse = "|")))
  )

#' Total number of files conforming
if (!exists("snakemake")) {
  sum(conformity$conforms)
  mean(conformity$conforms)
}

#' Total number of GEOs conforming
if (!exists("snakemake")) {
  conformity %>% 
    filter(conforms) %>% 
    pull(Accession) %>% 
    n_distinct()
}

#' Number of GEOs conforming per year
conformity_acc <- conformity %>% 
  group_by(Accession, year) %>% 
  summarise(
    conforms = case_when(
      any(conforms) ~ 1,
      TRUE ~ 0
    )) %>% 
  ungroup()
write_csv(conformity_acc, here("results/conformity_acc.csv"))

#' Total number of conforming GEO accessions.
#+
total_conforming <- conformity_acc %>% 
  count(conforms)

#' Conforming GEO submissions per year.
#+
if (!exists("snakemake")) {
  conformity_acc %>% 
    group_by(year) %>% 
    summarise(conforms = sum(conforms),
              n = n(),
              perc = conforms / n)
}

#' Number of sets with p-values, 
parsed_suppfiles_raw <- read_csv(here("results/data/parsed_suppfiles.csv"))
parsed_suppfiles <- parsed_suppfiles_raw %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

blacklist <- read_lines(here("results/data/blacklist.txt")) %>% 
  str_trim() %>% 
  tibble(suppfilenames = .) %>% 
  mutate(Accession = str_to_upper(str_extract(suppfilenames, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

c(parsed_suppfiles$Accession, blacklist$Accession) %>% 
  n_distinct()

unnested_suppfiles <- parsed_suppfiles %>% 
  mutate(
    splits = str_split(id, " from "),
    suppfile = map_chr(splits, ~ifelse(str_detect(.x[length(.x)], "RAW.tar$"), .x[length(.x) - 1], .x[length(.x)]))
  )

#' Number of p value sets per GEO submission
unnested_suppfiles %>% 
  filter(Type == "raw") %>% 
  count(Accession) %>% 
  summarise_at("n", list(mean = mean, median = median, min = min, max = max))

unnested_suppfiles %>% 
  filter(Type == "raw") %>% 
  count(Accession) %>% 
  summarise(n1 = mean(n == 1),
            n1_3 = mean(n <=3 ))

unnested_suppfiles_conforms <- unnested_suppfiles %>% 
  select(Accession, suppfile, note) %>% 
  distinct() %>% 
  mutate(
    conforms = !(is.na(suppfile) | str_detect(str_to_lower(suppfile), c("readme", "\\.[bs]am", "\\.bed") %>% str_c("(\\.gz)?", collapse = "|")))
  )

unnested_suppfiles_conforms %>% 
  group_by(Accession) %>% 
  summarise(
    conforms = case_when(
      any(conforms) ~ 1,
      TRUE ~ 0
    )) %>% 
  ungroup() %>% 
  summarise_at("conforms", mean)


#' Parse analysis platform
get_var <- function(x) {
  NAs <- is.na(x)
  if (all(NAs)) {
    return("unknown")
  }
  return(names(x)[!NAs])
}

get_n_tests <- function(x) {
  x %>% 
    str_remove_all("[\\[\\]]") %>% 
    str_split(", ", simplify = TRUE) %>% 
    as.numeric() %>% 
    sum(na.rm = TRUE)
}

pi0 <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  mutate(
    n_tests = map_dbl(hist, get_n_tests),
    prop_FDR_pval = ifelse(!is.na(pi0), FDR_pval / n_tests, NA_real_)
  ) %>% 
  select(Accession, id, Set, pi0, prop_FDR_pval) %>% 
  na.omit()

hist <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  select(Accession, id, Set, FDR_pval, hist)

# ##### 
# 
# 
# parsed_suppfiles_var <- parsed_suppfiles %>% 
#   left_join(acc_year) %>% 
#   select(Accession, id, Set, Type, Class, PDAT) %>% 
#   filter(!is.na(Type)) %>% 
#   pivot_wider(names_from = Type, values_from = Class) %>% 
#   group_by(id, Set) %>% 
#   mutate(var = get_var(c("logcpm"=logcpm, "rpkm"=rpkm, "fpkm"=fpkm, "basemean"=basemean, "aveexpr"=aveexpr))) %>%
#   ungroup() %>% 
#   mutate_at("Set", str_remove_all, '"') %>% 
#   rowwise()
# 
# platform_regex <- "deseq2?|de(g|x)seq|rockhopper|cuff(diff|links)|edger|clc(bio)? ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth|partek|(nrsa|nascent rna seq)|median ratio norm|rmats|ballgown|biojupie|seurat|exdega"
# analysis_platform_from_varset <- parsed_suppfiles_var %>% 
#   mutate(
#     analysis_platform_from_set = case_when(
#       str_detect(str_to_lower(Set), platform_regex) ~ str_extract(str_to_lower(Set), platform_regex), 
#       str_detect(str_to_lower(Set), "tagwise") ~ "edger.tagwise",
#       str_detect(str_to_lower(Set), "common") ~ "edger.common",
#       str_detect(str_to_lower(Set), "baggerley") ~ "clc genomics",
#       TRUE ~ NA_character_
#     ),
#     analysis_platform_from_var = case_when(
#       var == "basemean" & str_detect(str_to_lower(Set), "pval(?!ue)") ~ "deseq",
#       var == "basemean" & str_detect(str_to_lower(Set), "pvalue") ~ "deseq2",
#       var == "logcpm" ~ "edger",
#       var == "aveexpr" & str_detect(str_to_lower(Set), "p\\.value") & PDAT > "2014-01-01" ~ "limma",
#       var == "fpkm" & str_detect(str_to_lower(Set), "p_value") ~ "cuffdiff",
#       TRUE ~ NA_character_
#     ),
#     analysis_platform_from_filename = str_extract(str_to_lower(id), "deseq2?|de(g|x)seq|rockhopper|cuff(diff|links)|edger|clc ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth")
#   ) 
# 
# 
# analysis_platform_from_summaries <- function(x) {
#   x %>% 
#     map(str_to_lower) %>% 
#     map(str_extract_all, "deseq2?|de(g|x)seq|rockhopper|cuff(diff|links)|edger|clc ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth") %>% 
#     map(unlist) %>% 
#     map(unique)
# }
# 
# reorder <- function(x) {
#   str_c(sort(unlist(str_split(x, ";"))), collapse = ";")
# }
# 
# document_summaries_analysis_platform <- document_summaries %>% 
#   mutate(
#     analysis_platform_from_sum = map(summary, analysis_platform_from_summaries)
#   ) %>% 
#   select(Accession, analysis_platform_from_sum) %>% 
#   unnest(analysis_platform_from_sum) %>% 
#   unnest(analysis_platform_from_sum) %>% 
#   mutate(analysis_platform_from_sum = 
#            case_when(
#              analysis_platform_from_sum == "voom" ~ "limma",
#              TRUE ~ analysis_platform_from_sum)) %>% 
#   group_by(Accession) %>% 
#   summarise(analysis_platform_from_sum = str_c(unique(analysis_platform_from_sum), collapse = ";")) %>% 
#   ungroup() %>% 
#   rowwise() %>% 
#   mutate(analysis_platform_from_sum = reorder(analysis_platform_from_sum)) %>% 
#   ungroup()
# 
# accession_pubmedids <- document_summaries %>% 
#   select(Accession, PubMedIds) %>% 
#   drop_na() %>% 
#   distinct() %>% 
#   mutate(PubMedIds = str_split(PubMedIds, ";")) %>% 
#   unnest(PubMedIds)
# 
# detools_from_pmc_raw <- read_csv(here("results/detools_from_pmc.csv"), col_types = "ccc")
# detools_from_pmc <- detools_from_pmc_raw %>% 
#   left_join(accession_pubmedids) %>% 
#   select(Accession, PubMedIds, analysis_platform_from_pmc = detool)  %>% 
#   rowwise() %>% 
#   mutate(analysis_platform_from_pmc = reorder(analysis_platform_from_pmc)) %>% 
#   ungroup()
# 
# detools_from_pmc %>% 
#   select(-PubMedIds) %>% 
#   unique() %>% 
#   group_by(Accession) %>% 
#   summarise(analysis_platform_from_pmc = reorder(analysis_platform_from_pmc))
# 
# analysis_platform_lookup_table <- read_csv(here("results/analysis_platform_lookup_table.csv"))
# # analysis_platform_lookup_table[is.na(analysis_platform_lookup_table)] <- "na"
# pvalues1 <- analysis_platform_from_varset %>%
#   left_join(document_summaries_analysis_platform) %>% 
#   left_join(detools_from_pmc %>% select(-PubMedIds) %>% unique()) %>% 
#   # mutate_at(vars(starts_with("analysis_platform_from_")), ~ifelse(is.na(.x), "na", .x)) %>%
#   left_join(analysis_platform_lookup_table) %>% 
#   select(Accession, id, Set, Class = raw, analysis_platform_manual) %>% 
#   left_join(acc_year) %>% 
#   left_join(pi0) %>% 
#   left_join(hist) %>% 
#   select(-PDAT) %>% 
#   mutate(anticons = case_when(
#     Class == "anti-conservative" ~ 1,
#     TRUE ~ 0
#   )) %>% 
#   ungroup() %>% 
#   distinct()
# 
# pvalues1 %>% 
#   filter(analysis_platform_manual == "deseq2.var") %>% 
#   select(Accession) %>% 
#   left_join(document_summaries) %>% 
#   select(Accession, PubMedIds) %>% 
#   drop_na() %>% 
#   distinct() %>% 
#   pull(PubMedIds) %>% 
#   str_c(collapse = ",")
# 
# pmcid <- read_csv(here("results/pmcid.csv"), col_types = "ccc") %>% 
#   select(PMID, PMCID, DOI)
# 
# 
# dois_for_scraping <- pvalues1 %>% 
#   filter(str_detect(analysis_platform_manual, "\\.var") | is.na(analysis_platform_manual)) %>% 
#   select(Accession) %>% 
#   distinct() %>% 
#   left_join(document_summaries %>% 
#               select(Accession, PMID = PubMedIds) %>% 
#               drop_na() %>% 
#               inner_join(pmcid)) %>% 
#   filter(!is.na(PMID), !is.na(PMCID), !is.na(DOI))
# 
# dois_for_scraping %>% 
#   write_csv(here("results/dois_for_scraping.csv"))
# 
# 
# get_platform_matrix <- function(Accession, regex) {
#   nnn <- str_replace(Accession, "\\d{3}$", "nnn")
#   url <- glue::glue("https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{Accession}/matrix/{Accession}_series_matrix.txt.gz")
#   sm <- read_lines(url)
#   context <- sm[str_detect(str_to_lower(sm), regex)] %>% 
#     unique() %>% 
#     str_split("\\t") %>% 
#     map(unique) %>% 
#     flatten() %>% 
#     unique() %>% 
#     str_c(collapse = ";") %>% 
#     str_remove_all('(!Sample_data_processing;|")')
#   detool <- str_extract_all(str_to_lower(context), regex)[[1]] %>% 
#     unique() %>% 
#     str_c(collapse = ";")
#   tibble(detool = detool, context = context)
# }
# 
# safe_get_platform_matrix <- safely(get_platform_matrix)
# matrix_files2 <- dois_for_scraping %>% 
#   mutate(
#     res = map(Accession, safe_get_platform_matrix, regex = platform_regex)
#   )
# 
# detool_matrix <- matrix_files %>% 
#   mutate(res1 = map(res, "result")) %>% 
#   unnest(res1) %>% 
#   select(-res)
# 
# 
# detool_matrix %>% 
#   mutate_at(c("detool", "context"), ~ifelse(.x == "", NA_character_, .x)) %>% 
#   write_csv(here("results/detools_from_series_matrix.csv"))
# 
# detool_matrix %>% 
#   mutate_at(c("detool", "context"), ~ifelse(.x == "", NA_character_, .x)) %>% 
#   filter(is.na(detool)) %>% 
#   write_csv(here("results/detools_still_missing.csv"))
# 
# detool_matrix2 <- detool_matrix %>% 
#   mutate_at(c("detool", "context"), ~ifelse(.x == "", NA_character_, .x)) %>% 
#   filter(is.na(detool)) %>% 
#   select(-detool, -context) %>% 
#   mutate(
#     res = map(Accession, safe_get_platform_matrix, regex = str_c(platform_regex, "|partek|(nrsa|nascent rna seq)|median ratio norm|rmats|ballgown|biojupie|seurat|exdega"))
#   )
# 
# detool_matrix21 <- detool_matrix2 %>% 
#   mutate(res1 = map(res, "result")) %>% 
#   unnest(res1) %>% 
#   select(-res) %>% 
#   mutate_at(c("detool", "context"), ~ifelse(.x == "", NA_character_, .x))
# 
# 
# updated_regex <- list.files("results", pattern = "detools_pmcid_4\\d.csv", full.names = TRUE) %>% 
#   read_csv(col_types = "ccc") %>% 
#   bind_rows() %>% count(detool, sort=TRUE)
#   full_join(detool_matrix21)
# 
# more_platforms_update_regex <- updated_regex %>% 
#   filter(!is.na(detool)) %>% 
#   select(PMCID, detool)
# 
# more_platforms_update_regex %>% 
#   write_csv(here("results/detools_from_series_matrix_and_pmc_updated_regex.csv"))
#   
# 
# pvalues1 %>% 
#   filter(analysis_platform_manual == "deseq2.var") %>% 
#   group_by(Accession) %>% 
#   slice_sample(n = 1) %>% 
#   ungroup() %>% 
#   slice_sample(n = 10)
# 
# 
# 
# pvalues1 %>% 
#   ggplot(aes(analysis_platform_manual, pi0)) +
#   geom_violin() +
#   geom_jitter(size = 1/5) +
#   coord_flip()
# 
# 
# pvalues1 %>% 
#   ggplot(aes(y = analysis_platform_manual, fill = factor(anticons))) +
#   geom_bar() +
#   theme(axis.title.y = element_blank(),
#         legend.position = "bottom")
# 
# 
# pvalues1 %>% 
#   count(analysis_platform_manual, anticons) %>% 
#   group_by(analysis_platform_manual) %>% 
#   mutate(p = n / sum(n)) %>% 
#   filter(as.logical(anticons)) %>% 
#   ggplot() +
#   geom_col(aes(y = fct_reorder(analysis_platform_manual, p), x = p)) +
#   theme(axis.title.y = element_blank(),
#         legend.position = "bottom")
# 
# ###########

#' 
#' Analysis platform
#'  
#+
parse_analysis_platform_output <- read_csv(here("results/parse_analysis_platform_output.csv"))
pvalues <- parsed_suppfiles %>% 
  left_join(acc_year) %>% 
  select(Accession, id, Set, Type, Class, PDAT) %>% 
  filter(!is.na(Type)) %>% 
  pivot_wider(names_from = Type, values_from = Class) %>% 
  group_by(id, Set) %>% 
  mutate(var = get_var(c("logcpm"=logcpm, "rpkm"=rpkm, "fpkm"=fpkm, "basemean"=basemean, "aveexpr"=aveexpr))) %>%
  ungroup() %>% 
  mutate_at("Set", str_remove_all, '"') %>% 
  left_join(parse_analysis_platform_output) %>% 
  select(Accession, id, Set, Class = raw, analysis_platform) %>% 
  left_join(acc_year) %>% 
  left_join(pi0) %>% 
  left_join(hist) %>% 
  select(-PDAT) %>% 
  mutate(anticons = case_when(
    Class == "anti-conservative" ~ 1,
    TRUE ~ 0
  ))

write_csv(pvalues, here("results/pvalues.csv"))

#' Number of unique GEO ids imported
geo_import <- parsed_suppfiles %>% 
  filter(is.na(Type) | str_detect(Type, "raw"))
if (!exists("snakemake")) {
  geo_import %>% 
    pull(Accession) %>% 
    n_distinct()
}


#' Number of unique files downloaded
if (!exists("snakemake")) {
  geo_import %>% 
    select(Accession, id) %>% 
    mutate(file = str_remove(id, "^.+ from ")) %>% 
    pull(file) %>% 
    n_distinct()
}

#' Number of unique files imported
if (!exists("snakemake")) {
  geo_import %>% 
    pull(id) %>% 
    n_distinct()

parse_from <- function(x) {
  if (str_detect(x, "^sheet")) {
    nosheet <- str_split(x, " from ") %>% 
      unlist() %>% 
      str_subset("^[^sheet]") %>% 
      str_c(collapse = " from ")
    return(nosheet)
  }
  return(x)
}

geo_import %>% 
  mutate(id1 = map_chr(id, parse_from)) %>% 
  select(id1) %>% 
  distinct()
  
}

#' notes
if (!exists("snakemake")) {
  geo_import %>% 
    count(note) %>% 
    pull(note)
}

#' Files which we were able to import.
if (!exists("snakemake")) {
  geo_import %>% 
    mutate(imported = case_when(
      str_detect(str_to_lower(note), "error|i\\/o operation|not determine delim|codec can't decode|empty table|missing optional dependency") ~ "fail",
      is.na(note) ~ "yes",
      TRUE ~ "yes"
    )) %>% 
    group_by(Accession, id) %>% 
    summarise(
      imported = case_when(
        any(imported == "yes") ~ 1,
        TRUE ~ 0
      )) %>% 
    ungroup() %>% 
    count(imported)
}

#' P values
#+

#' Number of distinct GEO-s with p values
if (!exists("snakemake")) {
  n_distinct(pvalues$Accession)
}

#' Summary stats of p values per GEO
if (!exists("snakemake")) {
  pvalues %>% 
    count(Accession) %>% 
    summarise_at("n", list(mean = mean, median = median, max = max))
}

#' GEOs with single or 1-3 sets
#+
if (!exists("snakemake")) {
  pvalues %>% 
    count(Accession) %>% 
    count(n) %>% 
    mutate(p = nn / sum(nn),
           ps = cumsum(p))
}

#' Sample 1 p value set per accession
set.seed(11)
pvalues_sample <- pvalues %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  ungroup()
write_csv(pvalues_sample, here("results/pvalues_sample.csv"))

#' Save table ids for figure
#+
pvalues_sample %>% 
  select(id, Set) %>% 
  write_csv(here("results/pvalues_acc.csv"))

if (!exists("snakemake")) {
  pvalues_sample %>% 
    count(Class) %>% 
    summarise(Class, 
              n,
              p = signif(n / sum(n), 2))
}

