library(tidyverse)
library(lubridate)
library(cowplot)
library(patchwork)
library(here)
library(brms)
library(tidybayes)

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
# pp_check(m) + scale_x_log10()
p_year <- plot(conditional_effects(m), plot = FALSE)$year
# p_sigma <- plot(conditional_effects(m, dpar = "sigma"), plot = FALSE)$year +
#   scale_x_continuous(breaks = seq(from = min(data$year), to = max(data$year), 2)) +
#   theme_cowplot()
py <- p_year + 
  geom_point(
    data = drop_na(data), aes(year, log_reads), 
    inherit.aes = FALSE,
    position = position_jitter(0.5), size = 1/3
    ) +
  scale_x_continuous(breaks = seq(from = min(data$year), to = max(data$year), 2)) +
  scale_y_log10() +
  labs(y = "Median number of\nreads per sample, log10") +
  theme_cowplot()
py$layers <- rev(py$layers)

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
  labs(x = "Median number of\nreads per sample, log10",
       y = "Proportion of anti-conservative\np value sets") +
  theme_cowplot()

fig.cap <- str_wrap("Figure: Modelling of sequencing depth. A. 
                    Temporal change in sequencing depth. 
                    Log-transformed median libray size was modeled against year of GEO submission, model formula: log_reads ~ year, Student's likelihood. 
                    N = 31497. Blue line denotes model best fit. Gray area denotes 95% credible interval. Black points denote original data. Download model: https://gin.g-node.org/tpall/geo-htseq-paper/raw/96b9ccb224f18e4d0f0f23dabcd940a2011172be/models/spots__year.rds.
                    B. Proportion of anti-conservative p value histograms versus sequencing depth, model formula anticons ~ log_reads, bernoulli likelihood. N = 2081. Blue line denotes model best fit. Gray area denotes 95% credible interval. Download model: https://gin.g-node.org/tpall/geo-htseq-paper/raw/96b9ccb224f18e4d0f0f23dabcd940a2011172be/models/anticons__year_logreads.rds.",
                    width=120)
p <- py + pac + plot_annotation(tag_levels = "A", caption = fig.cap) 
ggsave(here("figures/figure_composition.png"), plot = p, width = 7, height = 5)
