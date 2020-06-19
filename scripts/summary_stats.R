library(tidyverse)
library(lubridate)

#' Number of unique GEOs
#+
document_summaries <- read_csv("data/document_summaries.csv")
document_summaries %>% 
  filter(PDAT<="2019-12-31") %>% 
  pull(Accession) %>% 
  n_distinct()

document_summaries %>% 
  filter(PDAT<="2019-12-31") %>% 
  group_by(year(PDAT)) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(cumsum = cumsum(n),
         perc = n / cumsum)

#' Supplementary files.
#+
suppfilenames <- read_lines("data/suppfilenames.txt") %>% 
  tibble(suppfilenames = .) %>% 
  filter(
    str_detect(suppfilenames, "suppl"),
    !str_detect(suppfilenames, "suppl/filelist.txt")
  )

acc_year <- document_summaries %>% 
  select(Accession, PDAT) %>% 
  mutate(year = year(PDAT))

conformity <- suppfilenames %>% 
  mutate(Accession = str_extract(suppfilenames, "GSE\\d+")) %>% 
  left_join(acc_year) %>% 
  mutate(
    conforms = !str_detect(str_to_lower(suppfilenames), 
                           "raw.tar|bam|sam|bed")
  ) 
#' Total number of GEOs conforming
conformity %>% 
  filter(conforms) %>% 
  pull(Accession) %>% 
  n_distinct()

#' Number of GEOs conforming per year
conformity %>% 
  group_by(year) %>% 
  count(conforms) %>% 
  mutate(
    cumsum = cumsum(n),
    perc = n / cumsum
  ) %>% 
  filter(conforms)

#' Number of sets with p-values
parsed_suppfiles <- read_csv("data/parsed_suppfiles.csv") %>% 
  filter(is.na(Type) | str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())
parsed_suppfiles %>% 
  pull(Accession) %>% 
  n_distinct()

parsed_suppfiles %>% 
  filter(str_detect(Type, "raw")) %>% 
  pull(Accession) %>% 
  n_distinct()
  
