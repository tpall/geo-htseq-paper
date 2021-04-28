if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}

#+ libs
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(lubridate)
library(cowplot)
library(Hmisc)
library(brms)
library(rstan)
library(here)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

#+ params
is_ci <- function() {
  "CI" %in% Sys.getenv()
}
chains <- ifelse(is_ci(), 1, 4)
cores <- chains
refresh = 0
rstan_options(auto_write = TRUE, javascript = FALSE)

#' Importing publication data.
#+
publications <- read_csv(
  here("resources/data/publications.csv"),
  col_types = cols(
    .default = col_character(),
    Id = col_character(),
    HasAbstract = col_double(),
    PmcRefCount = col_double()
  ))
pubs <- publications %>% 
  mutate_at(c("ISSN", "ESSN"), str_remove_all, "-") %>% 
  mutate(year = as.numeric(str_extract(PubDate, "^\\d{4}"))) %>% 
  select(PubMedIds = Id, year, title = FullJournalName, ISSN, ESSN) %>% 
  pivot_longer(cols = c("ISSN", "ESSN")) %>% 
  select(-name) %>% 
  drop_na()

#' Importing citaion data.
#+
citations <- read_csv(
  here("resources/data/scopus_citedbycount.csv"),
  col_types = cols(
    PubMedIds = col_character(),
    citations = col_double()
  ))

#' Import Citescore data.
#+
path <- here("resources/data/CiteScore 2011-2019 new methodology - October 2020.xlsx")
citescore <- excel_sheets(path) %>% 
  str_subset("^CiteScore")
citescores_raw <- citescore %>% 
  map(~read_xlsx(path, sheet = .x)) %>% 
  set_names(citescore) %>% 
  bind_rows(.id = "year")
citescores <- citescores_raw %>% 
  select(year, Title, CiteScore, ISSN = `Print ISSN`, ESSN = `E-ISSN`) %>% 
  mutate_at("year", str_remove, "CiteScore ") %>% 
  mutate_at("year", as.numeric) %>% 
  pivot_longer(cols = c("ISSN", "ESSN")) %>% 
  select(-name) %>% 
  drop_na()

#' Merge CiteScores with publication data.
#+
pubs_all <- pubs %>% 
  left_join(citescores) %>% 
  select(PubMedIds, year, Title, CiteScore) %>% 
  distinct() %>% 
  filter(year <= 2019)

pubs_citescore <- pubs_all %>% 
  drop_na() %>% 
  left_join(citations)

#' Importing document summaries.
#+
document_summaries <- read_csv(
  here("resources/data/document_summaries.csv"),
  col_types = cols(
    .default = col_character(),
    Id = col_character(),
    GDS = col_logical(),
    GSE = col_double(),
    ptechType = col_logical(),
    valType = col_logical(),
    SSInfo = col_logical(),
    subsetInfo = col_logical(),
    PDAT = col_date(format = ""),
    n_samples = col_double(),
    SeriesTitle = col_logical(),
    PlatformTitle = col_logical(),
    PlatformTaxa = col_logical(),
    SamplesTaxa = col_logical()
  ))

doc_pubs <- document_summaries %>% 
  select(Accession, PDAT, PubMedIds) %>% 
  mutate(published = as.numeric(!is.na(PubMedIds)),
         year = year(PDAT)) %>% 
  group_by(year) %>% 
  summarise(published = sum(published), total = n())

cbind(year = doc_pubs$year, binconf(doc_pubs$published, doc_pubs$total)) %>% 
  as_tibble() %>% 
  ggplot() +
  geom_line(aes(year, PointEst)) +
  geom_ribbon(aes(year, ymin = Lower, ymax = Upper), alpha = 0.3) +
  scale_x_continuous(breaks = doc_pubs$year) +
  labs(
    x = "Year", 
    y = "Proportion of GEO series with journal article",
    caption = "Proportion of published GEO series that are related to journal publication. Gray area denotes binomial 95% confidence interval.")

#' Add Accession to PubMedIds
#+
acc_pmid <- document_summaries %>% 
  select(Accession, PubMedIds, taxon) %>% 
  mutate(PubMedIds = map(PubMedIds, str_split, ";"),
         PubMedIds = map(PubMedIds, 1)) %>% 
  unnest(PubMedIds) %>% 
  drop_na()

pvalues_sample <- read_csv(here("results/pvalues_sample.csv"))
data <- pvalues_sample %>% 
  select(Accession, Class) %>% 
  left_join(acc_pmid %>% 
              left_join(pubs_citescore) %>% 
              drop_na()) %>% 
  drop_na() %>% 
  mutate(anticons = as.numeric(Class == "anti-conservative"))

taxa <- data %>% 
  count(taxon) %>% 
  arrange(desc(n)) %>% 
  head(20) %>% 
  pull(taxon)

data <- data %>% 
  mutate(Taxon = case_when(
    taxon %in% taxa ~ taxon,
    TRUE ~ "other"
  ))

f <- anticons ~ CiteScore + Taxon + year
family <- bernoulli()
mod <- brm(
  formula = f, 
  data = data,
  family = family, 
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = here("results/models/anticons__CiteScore_Taxon_year.rds")
  )
pp_check(mod)
conditions <- make_conditions(data, "Taxon")
p <- plot(
  conditional_effects(
    mod, 
    effects = "CiteScore", 
    conditions = conditions), 
  line_args = list(color = "black"), plot = FALSE)
p$CiteScore + facet_wrap(~Taxon)
