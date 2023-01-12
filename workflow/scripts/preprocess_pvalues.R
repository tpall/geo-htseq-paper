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


#' Number of sets with p-values, 
parsed_suppfiles_raw <- read_csv(here("results/data/parsed_suppfiles.csv"))
parsed_suppfiles <- parsed_suppfiles_raw %>% 
  mutate(
    Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  mutate_at("Set", str_remove_all, '"') %>% 
  select(Accession, everything())

#' Preparsed analysis platforms
parse_analysis_platform_output <- read_csv(here("results/parse_analysis_platform_output.csv"))

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
    str_split(",", simplify = TRUE) %>% 
    str_trim() %>% 
    as.numeric() %>% 
    sum(na.rm = TRUE)
}

pi0 <- parsed_suppfiles %>% 
  filter(Type == "raw", !is.na(Type)) %>% 
  mutate(
    n_tests = map_dbl(hist, get_n_tests),
    prop_FDR_pval = ifelse(!is.na(pi0), FDR_pval / n_tests, NA_real_)
  ) %>% 
  select(Accession, id, Set, pi0, prop_FDR_pval) 

hist <- parsed_suppfiles %>% 
  filter(Type == "raw", !is.na(Type)) %>% 
  select(Accession, id, Set, FDR_pval, hist)


#' 
#' Analysis platform
#'  
#+
pvalues <- parsed_suppfiles %>% 
  left_join(acc_year) %>% 
  select(Accession, id, Set, Type, Class, PDAT) %>% 
  filter(!is.na(Type)) %>% 
  pivot_wider(names_from = Type, values_from = Class) %>% 
  group_by(id, Set) %>% 
  mutate(var = get_var(c("logcpm"=logcpm, "rpkm"=rpkm, "fpkm"=fpkm, "basemean"=basemean, "aveexpr"=aveexpr))) %>%
  ungroup() %>% 
  # mutate_at("Set", str_remove_all, '"') %>% 
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
pvalues %>% 
  filter(as.logical(anticons)) %>% 
  pull(id) %>% 
  str_split(" from ") %>% 
  map(tail, 1) %>% 
  write_lines(here("results/anti_conservative_suppfiles.csv"))

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


if (!exists("snakemake")) {
  pvalues_sample %>% 
    count(Class) %>% 
    summarise(Class, 
              n,
              p = signif(n / sum(n), 2))
}

