library(tidyverse)
library(here)
sc_acc <- read_lines(here("resources/data/single-cell.csv"))
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample %>% 
  filter(Accession %in% sc_acc) %>% 
  ggplot() +
  geom_bar(aes(de_tool))
