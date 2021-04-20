library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(lubridate)
library(here)
library(networkD3)

#' Number of sets with p-values
#+
col_types <- cols(
  id = col_character(),
  Type = col_character(),
  Class = col_character(),
  Conversion = col_character(),
  pi0 = col_double(),
  FDR_pval = col_double(),
  hist = col_character(),
  note = col_character(),
  Set = col_character()
)
parsed_suppfiles_raw <- read_csv(here("results/data/parsed_suppfiles.csv"), col_types = col_types)
parsed_suppfiles <- parsed_suppfiles_raw %>% 
  filter(!str_detect(id, "_RAW.tar")) %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

set.seed(11)
pvalues_sample <- parsed_suppfiles %>% 
  filter(!is.na(Conversion)) %>% 
  mutate(Type = if_else(Type == "raw", "raw", "filtered")) %>% 
  select(Accession, id, Type, Class, Set) %>% 
  pivot_wider(names_from = Type, values_from = Class) %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  ungroup()

links <-  pvalues_sample %>% 
  count(raw, filtered, name = "value") %>% 
  mutate(source = str_c(raw, "_raw"),
         target = str_c(filtered, "_filtered")) %>% 
  select(source, target, value)
nodes <- data.frame(name = unique(c(links$source, links$target)), 
                    stringsAsFactors = FALSE)

# change the source and target variables to be the zero-indexed position of
# each node in the new nodes data frame
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1
nodes <- nodes %>% 
  mutate(name = str_remove(name, "_.*$"))

sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", 
              fontSize = 24, nodeWidth = 30, fontFamily = "Helvetica",
              height = 600, width = 750)
