library(tidyverse)
library(lubridate)
library(cowplot)
library(patchwork)
library(here)
library(brms)

old <- theme(theme_cowplot())

col_types <- cols(
  geo_accession = col_character(),
  study_accession = col_character(),
  run_accession = col_character(),
  sample_name = col_character(),
  sample_title = col_character(),
  spots = col_double(),
  bases = col_double(),
  tax_id = col_character(),
  organism = col_character(),
  LIBRARY_STRATEGY = col_character(),
  LIBRARY_SOURCE = col_character(),
  LIBRARY_SELECTION = col_character(),
  LIBRARY_LAYOUT = col_character(),
  PLATFORM = col_character(),
  MODEL = col_character(),
  err = col_character()
)
spots_raw <- read_csv(here("resources/data/spots.csv"), col_types = col_types)
spots <- spots_raw %>% 
  group_by(Accession = geo_accession) %>% 
  summarise_at("spots", median)
col_types <- cols(
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
)
doc_sums <- read_csv(here("resources/data/document_summaries.csv"), col_types = col_types)
composition <- doc_sums %>% 
  select(Accession, PDAT) %>% 
  left_join(spots) %>% 
  drop_na() %>% 
  mutate(year = year(PDAT))

composition %>% 
  ggplot() +
  geom_histogram(aes(spots)) +
  scale_x_log10()

data <- composition %>% 
  mutate(log_reads = log10(spots))
  
f <- log_reads ~ year
m <- brm(f, 
         data = data, 
         family = student(), 
         file = here("results/models/spots__year.rds")
)

p <- plot(conditional_effects(m), plot = FALSE, points = TRUE)
py <- p$year + 
  scale_x_continuous(breaks = seq(from = min(data$year), to = max(data$year), 2)) +
  theme_cowplot()

col_types <- cols(
  Accession = col_character(),
  id = col_character(),
  Set = col_character(),
  Class = col_character(),
  analysis_platform = col_character(),
  year = col_double(),
  pi0 = col_double(),
  FDR_pval = col_double(),
  hist = col_character(),
  anticons = col_double()
)

pvalues_sample <- read_csv(here("results/pvalues_sample.csv"), col_types = col_types)
composition %>% 
  mutate(log_reads = log10(spots))
data <- pvalues_sample %>% 
  left_join(data)
f <- anticons ~ year + log_reads
m <- brm(f, 
         data = data, 
         family = bernoulli(), 
         file = here("results/models/anticons__year_logreads.rds")
)

p <- plot(conditional_effects(m), plot = FALSE)
pac <- p$log_reads +
  theme_cowplot()

p <- py + pac + plot_annotation(tag_levels = "A")
ggsave(here("figures/figure_composition.png"), plot = p)

