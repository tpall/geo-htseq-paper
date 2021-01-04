
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
  
}

#' Import required libraries.
#+ libs
library("readr")
library("stringr")
library("tidyr")
library("here")
library("dplyr")
library("purrr")
library("skimr")

#' Helper function to write output description to README.
#+
write_desc <- function(object, file, file_desc, var_desc, append_to_file = NULL) {
  skimmed <- skim(object)
  skim_var <- skimmed$skim_variable
  skim_type <- skimmed$skim_type
  variables <- paste(skim_var, skim_type, sep=": ") %>% 
    paste(var_desc, sep = ", ") %>% 
    paste(collapse = "\n\t\t") %>% 
    paste0("\t\t", .) %>% 
    paste0("\tvariables:\n", .)
  desc <- paste(basename(file), file_desc, sep = ": ") %>% 
    paste(variables, sep = "\n")
  write(paste0("\n",desc), file = here(append_to_file), append = TRUE)
}

#' Document summaries.
#+
document_summaries <- read_csv(here("results/data/document_summaries.csv"), col_types = cols(Id = col_character()))
docsums <- document_summaries %>% 
  select(Accession, title, summary, taxon, PDAT, n_samples, PubMedIds)

#' Supplementary file names.
#+
suppfile_names <- read_delim(here("results/data/suppfilenames.txt"), delim = " ", col_names = "suppfile") %>% 
  mutate_at("suppfile", str_remove, "suppl/") %>% 
  mutate(Accession = str_extract(suppfile, "GSE\\d+")) %>% 
  fill(Accession) %>% 
  group_by(Accession) %>% 
  summarise(n_suppfiles = n(), 
            suppfiles = str_c(suppfile, collapse = "; "))

#' Merge supplementary file names with document summaries
#+
geo_suppfiles <- docsums %>% 
  left_join(suppfile_names) %>% 
  select(Accession, PDAT, title, summary, PubMedIds, n_samples, taxon, n_suppfiles, suppfiles)

#' Write merged file to csv.
#+
geo_suppfiles %>% 
  write_csv(here("results/geo_suppfiles.csv"))

#' Import parsed supplementary files.
#+
parsed_suppfiles <- read_csv(here("results/data/parsed_suppfiles.csv"))
suppfiles <- parsed_suppfiles %>% 
  filter(is.na(Type) | Type == "raw") %>% 
  select(id, Class, pi0, FDR_pval, hist, note, Set) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())

#' Print import messages. NA -- p values present.
#+
suppfiles %>% 
  count(note) %>% 
  arrange(desc(n))

#' Import publications associated with GEO submissions.
#+
publications <- read_csv("results/data/publications.csv", col_types = cols(Id = col_character()))
pubs <- publications %>% 
  select(PubMedIds = Id, PubDate, EPubDate, Source, AuthorList, LastAuthor, Title, DOI)

#' Import citation data.
#+
scopus_citedbycount <- read_csv("results/data/scopus_citedbycount.csv", col_types = cols(PubMedIds = col_character()))
pubs <- pubs %>% 
  left_join(scopus_citedbycount)

pubmedids<- geo_suppfiles %>% 
  select(Accession, PubMedIds) %>% 
  rowwise() %>% 
  mutate(PubMedIds = map(PubMedIds, str_split, pattern = ";")) %>% 
  unnest(PubMedIds) %>% 
  unnest(PubMedIds)

geo_publications <- pubmedids %>% 
  left_join(pubs)

#' Merge imported supplementary files data with p values.
#+
pvalues <- read_csv(here("results/pvalues.csv"))
imported_suppfiles <- suppfiles %>% 
  left_join(select(pvalues, id, Set, analysis_platform)) %>% 
  mutate(note = case_when(
    is.na(note) ~ "OK",
    TRUE ~ note
  ))
imported_suppfiles %>% 
  write_csv(here("results/imported_suppfiles.csv"))

#' Write dataset descriptions to README
#+
file.remove(here("results/README.txt"))
var_desc <- c(
  "GEO accession",
  "GEO title",
  "GEO summary",
  "PubMed id number(s)",
  "taxon(s) used in a study",
  "supplementary file names",
  "publication date",
  "number of samples",
  "number of supplementary files"
)
write_desc(geo_suppfiles, 
           "results/geo_suppfiles.csv", 
           "document summary and supplementary files of GEO submissions", 
           var_desc, 
           "results/README.txt")

var_desc <- c(
  "GEO accession",
  "supplementary table id, generated from file name, sheet name (optional), archive name (optional)",
  "class assigned to p value set",
  "list, number of p values belonging to equally spaced bins between 0,1",
  "supplementary file import message",
  "p value column name identified from supplementary table",
  "differential expression analysis platform inferred from column name pattern",
  "proportion of true null hypotheses estimated from p value set",
  "number of p values below FDR threshold (<0.05)"
)
write_desc(imported_suppfiles, 
           "results/imported_suppfiles.csv", 
           "import info and raw p value sets identified from downloaded GEO supplementary files", 
           var_desc, 
           "results/README.txt")

var_desc <- c(
  "GEO accession",
  "Pubmed id number",
  "publication date",
  "electronic publication date",
  "publication source title",
  "authors list",
  "last author name",
  "article title",
  "Digital Object Identifier",
  "number of citations"
)
write_desc(geo_publications, 
           "results/geo_publications.csv", 
           "publications citing GEO submissions", 
           var_desc, 
           "results/README.txt")

sequencing_metadata <- read_csv(here("results/sequencing_metadata.csv"))
var_desc <- c(
  "GEO accession",
  "sequencing library type",
  "sequencing library source",
  "approach for sequencing library selection",
  "sequencing library layout",
  "sequencing instrument platform",
  "sequencing instrument model",
  "source taxon id",
  "number of reads in library"
)
write_desc(sequencing_metadata, 
           "results/sequencing_metadata.csv", 
           "sequencing metadata from SRA", 
           var_desc, 
           "results/README.txt")
