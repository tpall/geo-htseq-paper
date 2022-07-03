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

pi0 <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  select(Accession, id, Set, pi0) %>% 
  na.omit()

hist <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  select(Accession, id, Set, FDR_pval, hist)

pvalues <- parsed_suppfiles %>% 
  left_join(acc_year) %>% 
  select(Accession, id, Set, Type, Class, PDAT) %>% 
  filter(!is.na(Type)) %>% 
  pivot_wider(names_from = Type, values_from = Class) %>% 
  group_by(id, Set) %>% 
  mutate(var = get_var(c("logcpm"=logcpm, "rpkm"=rpkm, "fpkm"=fpkm, "basemean"=basemean, "aveexpr"=aveexpr))) %>%
  ungroup() %>% 
  mutate( # parsing analysis platform using expression level variable name
    analysis_platform_from_expression = case_when(
      var == "basemean" ~ "deseq",
      (var == "aveexpr" & PDAT > "2014-01-01") ~ "limma",
      var == "logcpm" ~ "edger",
      var == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
      TRUE ~ "unknown"
    ), # parsing analysis platform using file names
    analysis_platform_from_filename = case_when(
      str_detect(str_to_lower(id), "deseq") ~ "deseq",
      str_detect(str_to_lower(id), "edger") ~ "edger",
      str_detect(str_to_lower(id), "limma") ~ "limma",
      str_detect(str_to_lower(id), "cuff") ~ "cuffdiff",
      TRUE ~ "unknown"
    ), # assigning analysis platform using file names and expression level variable name
    analysis_platform = case_when(
      analysis_platform_from_expression == analysis_platform_from_filename ~ analysis_platform_from_expression,
      analysis_platform_from_filename == "unknown" ~ analysis_platform_from_expression,
      TRUE ~ analysis_platform_from_filename
      )
  ) %>%
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

