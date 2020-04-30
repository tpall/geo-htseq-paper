library("tidyverse")
library("sparkline")
library("formattable")
library("here")

document_summaries <- read_csv("data/document_summaries.csv") %>% 
  mutate(sample_acc = str_extract_all(Samples, "GSM\\d+")) %>% 
  select(accession = Accession, sample_acc) %>% 
  unnest(cols = sample_acc)
pvals_wide <- read_csv("output/pvals_wide.csv")
spots <- read_csv("data/spots.csv") %>% 
  mutate(sample_acc = str_replace(name, "(^GSM\\d+).+", "\\1"),
         name = str_replace(name, "(^GSM\\d+: )(.+)", "\\2")) %>% 
  separate(name, c("title", "species", "exp_type"), sep = "; ") %>% 
  select(sample_acc, starts_with("total"))
spot_merged <- inner_join(spots, document_summaries)
spot_summary <- spot_merged %>% 
  group_by(acc) %>% 
  summarise_at(vars(starts_with("total")), list(mean = mean, median = median, sd = sd, min = min, max = max, n = length))
