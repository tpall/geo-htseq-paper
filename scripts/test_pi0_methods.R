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

supp_smo <- read_csv("data/parsed_suppfiles.csv")
supp_lfdr <- read_csv("~/Downloads/parsed_suppfiles.csv")
document_summaries <- read_csv("data/document_summaries.csv")

gses <- setdiff(supp_smo$id, supp_lfdr$id) %>% 
  str_extract("GSE\\d+") %>% 
  unique()
length(gses)

document_summaries %>% 
  filter(Accession %in% gses) %>% 
  pull(PDAT) %>% 
  range()

supp_smo %>% 
  filter(Type == "raw", !is.na(Set))
supp_lfdr %>% 
  filter(Type == "raw", !is.na(Set))

anti_conservative <- supp_smo %>% 
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

anti_conservative_lfdr <- supp_lfdr %>% 
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

