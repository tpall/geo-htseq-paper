library("tidyverse")
library("sparkline")
library("formattable")
library("here")

#' ## P value stats
#' Import dataset
imported <- read_csv(here("output/parsed_suppfiles.csv"))
imported <- imported %>% 
  mutate(accession = str_extract(id, "GSE\\d+")) %>% 
  select(accession, everything())

#' Pivot wide
pvalues <- imported %>% 
  filter(!is.na(Type)) %>% 
  mutate(Class = str_replace(Class, "conservative", "cons."),
         Conversion = str_replace(Conversion, "improvement", "impr."))
before <- pvalues %>% filter(Type == "raw")
after <- pvalues %>% filter(Type != "raw") %>% 
  rename_at(vars(Class, pi0, FDR_pval, hist), str_c, "_after") %>% 
  rename(filter_var = Type)
wide <- full_join(before, after) %>% 
  select(accession, id, starts_with("class"), Conversion, starts_with("pi0"), starts_with("hist"), filter_var, Set)


#' Generate html table: beware of very large table!
# tab <- wide %>% 
#   mutate_at(vars(starts_with("hist")), ~map(.x, ~as.integer(unlist(str_extract_all(.x, "\\d+"))))) %>% 
#   mutate_at(vars(starts_with("hist")), ~map(.x, spk_chr, type = "bar")) %>%
#   formattable() %>%
#   formattable::as.htmlwidget() %>%
#   spk_add_deps()

#' Conversion
#+ echo=FALSE
wide %>% 
  count(Class, Class_after) %>% 
  group_by(Class) %>% 
  mutate(`%` = formatC((n / sum(n[!is.na(Class_after)]) * 100), digits = 1, format = "f"),
         `%` = ifelse(is.na(Class_after), NA, `%`)) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

#' Number of p-value sets per GEO accession
#+
wide %>% 
  count(accession) %>% 
  arrange(desc(n))


#' ## Experiment metadata
#' Importing experiment metadata (dates and etc) 
document_summaries <- read_csv("output/document_summaries.csv")
range(document_summaries$PDAT)

#' Number of unique GEO sets during period.
#+
n_distinct(document_summaries$Accession)

#' Number of GEO sets with supplementary files that were imported.
#+
n_distinct(imported$accession)

#' Number of files with GEO sets.
#'+
imported %>% 
  filter(note == "no pvalues" | Type == "raw") %>% 
  select(accession, id, Type, note)


#' Importing publications metadata
publications <- read_csv("output/publications.csv", col_types = str_c(rep("c", 26), collapse = "")) %>% 
  rename(PubMedIds = Id)

#' Single-cell experiments
single_cell <- read_csv("output/single-cell.csv")

#' Citations
citations <- read_csv("output/scopus_citedbycount.csv")

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
  select(accession = Accession, PDAT, taxon) %>% 
  inner_join(wide) %>% 
  mutate(single_cell = accession %in% single_cell$Accession) %>% 
  mutate(platform = case_when(
    filter_var == "basemean" ~ "deseq",
    filter_var == "aveexpr" ~ "limma",
    filter_var == "logcpm" ~ "edger",
    filter_var == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
    TRUE ~ "unknown"
  ))
range(pvals_wide$PDAT)

pvals_wide %>% 
  count(platform)

pvals_wide %>% 
  count(taxon) %>% 
  arrange(desc(n))
