library(tidyverse)
library(data.table)
library(limma)

files <- list.files("~/Downloads/output_12/suppl/", full.names = TRUE, recursive = TRUE)
safe_fread <- safely(fread)
imported <- files %>% 
  map(fread)

tibbles <- imported %>% 
  map(as_tibble, .name_repair = "universal")
  
pvalues <- tibbles %>% 
  map(select, matches("p[^a-zA-Z]{0,4}val")) %>% 
  map(select, -matches("adj|fdr|corr|thresh"))

names(pvalues) <- basename(files)
pi0 <- pvalues %>% 
  map(summarise_all, propTrueNull, method = "hist")

anti_conservative <- read_csv("data/parsed_suppfiles.csv") %>% 
  filter(Type == "raw", !is.na(pi0)) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())

files_to_check <- anti_conservative %>% 
  filter(near(pi0, 1)) %>% 
  select(id, pi0, FDR_pval, Set)

anti_conservative %>% 
  ggplot() +
  geom_histogram(aes(pi0), binwidth = 0.05) +
  labs(title = "smoother")

anti_conservative_lfdr <- read_csv("~/Downloads/parsed_suppfiles.csv") %>% 
  filter(Type == "raw", !is.na(pi0)) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())

anti_conservative_lfdr %>% 
  ggplot() +
  geom_histogram(aes(pi0), binwidth = 0.05) +
  labs(title = "lfdr")

list(anti_conservative, anti_conservative_lfdr) %>% 
  set_names(c("smoother", "lfdr")) %>% 
  bind_rows(.id = "method") %>% 
  ggplot() +
  geom_density(aes(pi0, linetype = method))

