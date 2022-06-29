
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}

#+ libs
library(stats) # masks filter
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(grid)
library(gridExtra)
library(gtable)
library(extrafont)
library(cowplot)
library(patchwork)
library(viridis)
library(brms)
library(rstan)
library(tidybayes)
library(modelr)
library(magick)
library(networkD3)
library(webshot)
library(here)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))


if (!is_phantomjs_installed()) {
  install_phantomjs()
} 

#+ params
is_ci <- function() {
  "CI" %in% Sys.getenv()
}
chains <- 1 #ifelse(is_ci(), 1, 4)
cores <- chains
refresh = 200
rstan_options(auto_write = TRUE, javascript = FALSE)
if (!dir.exists("results/models")) {
  dir.create("results/models", recursive = TRUE)
}

#+ data
conformity_acc <- read_csv(here("results/conformity_acc.csv"))
pvalues <- read_csv(here("results/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_filtered <- read_csv(here("results/pvalues_filtered.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_filtered_sample <- read_csv(here("results/pvalues_filtered_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))

#'
#+ fig1
f <- conforms ~ year
# We will rescale years, so that year = 0 is the year of first submission
data <- conformity_acc %>% 
  mutate_at("year", ~.x - min(.x))
family <- bernoulli()
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = prior("normal(0, 2)", class = "b"),
           refresh = refresh,
           chains = chains,
           cores = cores,
           control = list(adapt_delta = 0.95, max_treedepth = 12),
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/conforms_year.rds"),
           file_refit = "on_change")

p <- plot(conditional_effects(mod), plot = FALSE, line_args = list(color = "black", size = 1/2))$year
p <- p + 
  geom_point(
    data = data %>% group_by(year) %>% summarise_at("conforms", list(conforms = mean)), 
    aes(year, conforms), 
    inherit.aes = FALSE,
    size = 1/2) +
  labs(x = "Year", y = "Proportion of submissions conforming\nwith GEO submission guidelines") +
  scale_x_continuous(breaks = seq(0, 12, by = 2), labels = seq(2008, 2020, by = 2))
ggsave(here("figures/figure_1.pdf"), plot = p, height = 6, width = 7, dpi = 300, units = "cm")

#' 
#+ fig2
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
  filter(str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())


pvalues_acc <- read_csv(here("results/pvalues_acc.csv"))
suppfiles_sample <- parsed_suppfiles %>% 
  inner_join(pvalues_acc)


files_to_check <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  filter(near(pi0, 1)) %>% 
  select(id, pi0, FDR_pval, Set) %>% 
  pull(id) %>% 
  str_split(" from ") %>% 
  flatten_chr()


files_to_check[str_detect(files_to_check, "^GSE")] %>% 
  unique() %>% 
  str_c("output/suppl/", .) %>% 
  write_lines(here("results/files_to_check_pi0.txt"))


qc_threshold <- function(x, fdr) {
  bins <- length(x)
  qbinom(1 - 1 / bins * fdr, sum(x), 1 / bins)
}

hist_examples <- list(
  Class = c("conservative", "anti-conservative", "uniform", "bimodal", "other"), 
  id = c("GSE63555_Gene_expression_cuffdiff_fpkm.txt.gz", "GSE89511_diff.geneFpkm.exon.glm.LogFc.0.exc.Gapdh.RPMI8226_MCL1.txt.gz", "GSE102826_model_fc.txt.gz", "GSE115649_P20WT_P20M_results.csv.gz", "GSE98869_diffexp.tsv.gz"),
  Set = c("p_value", "vec.v.cd86.pvalue", "p.value (mutunc5 - control)", "pvalue", "pvalue")
) %>% 
  as_tibble()

hist_data <- parsed_suppfiles %>% 
  inner_join(hist_examples)

hist_data_plots <- hist_data %>% 
  mutate(hist = str_remove_all(hist, "[:punct:]"),
         hist = str_split(hist, " "),
         hist = map(hist, as.numeric)) %>% 
  select(Accession, id, Class, hist) %>%
  mutate(QC_thr = map_dbl(hist, qc_threshold, fdr = 0.05))


class_counts <- suppfiles_sample %>% 
  count(Class, name = "N")

fit <- brm(Class ~ 1, 
           data = suppfiles_sample, 
           family = categorical(), 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/Class_1.rds"),
           file_refit = "on_change"
)

pe <- posterior_epred(fit)
classes_props <- pe[1:dim(pe)[1], 1, 1:5] %>% 
  as_tibble() %>% 
  map(mean_hdi) %>% 
  bind_rows(.id = "Class") %>% 
  mutate_at(vars(y, ymin, ymax), signif, digits = 2) %>% 
  mutate(`Fraction [95% CI]` = str_c(y, " [", ymin, "; ", ymax, "]"),
         Example = NA) %>% 
  select(Class, `Fraction [95% CI]`, Example)

plot_data <- hist_data_plots %>% 
  select(Class) %>% 
  left_join(class_counts) %>% 
  left_join(
    classes_props
  ) %>% 
  arrange(desc(N)) %>% 
  select(Class, N, `Fraction [95% CI]`) %>% 
  ungroup()

fig2a <- tableGrob(
  arrange(plot_data, Class), 
  rows = NULL, 
  theme = ttheme_minimal(base_size = 8))
fig2a <- gtable_add_grob(fig2a,
                         grobs = segmentsGrob( # line across the bottom
                           x0 = unit(0,"npc"),
                           y0 = unit(0,"npc"),
                           x1 = unit(1,"npc"),
                           y1 = unit(0,"npc"),
                           gp = gpar(lwd = 2.0)),
                         t = 1, b = 1, l = 1, r = 3)

fig2b <- hist_data_plots %>% 
  mutate(
    data = map(hist, ~tibble(x = seq(0, 1, length.out = length(.x)), y = .x)),
    Class = factor(Class, levels = c("uniform", "bimodal", "other", "anti-conservative", "conservative"))
  ) %>% 
  select(Class, QC_thr, data) %>% 
  unnest(data) %>% 
  ggplot() +
  geom_col(aes(x, y)) +
  geom_hline(aes(yintercept = QC_thr), hist_data_plots, color = "red") +
  facet_wrap(~Class, ncol = 1,scales = "free") +
  theme_minimal_vgrid(font_size = 8) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

layout <- "
A###
ABBB
A###
"
p2 <- wrap_elements(fig2b)  + wrap_elements(fig2a) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout)
ggsave(here("figures/figure_2.pdf"), plot = p2, width = 12, height = 8, units = "cm", dpi = 300)
ggsave(here("figures/figure_2.eps"), plot = p2, width = 12, height = 8, units = "cm", dpi = 300)
ggsave(here("figures/figure_2.png"), plot = p2, width = 12, height = 8, units = "cm", dpi = 300)

#'
#+ fig3
f <- Class ~ year + (year | de_tool)
family <- categorical()
# We will rescale years, so that year = 0 is the year of first submission
data <- pvalues_sample %>% 
  mutate_at("year", ~.x - min(.x)) %>%
  select(Class, year, de_tool) %>% 
  drop_na()
get_prior(f, data, family)

data %>% 
  count(de_tool, Class, year) %>% 
  group_by(de_tool, year) %>% 
  mutate(p = n / sum(n)) %>% 
  ggplot(aes(year, p, color = Class)) +
  geom_line() +
  geom_point() +
  facet_wrap(~de_tool)

priors <- c(
  set_prior("lkj(3)", class = "cor"),
  set_prior("student_t(3, 0, 1)", class = "Intercept"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "muuniform"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "muuniform"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "mubimodal"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "mubimodal"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "muconservative"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "muconservative"),
  set_prior("normal(0, 0.1)", class = "b", dpar = "muother"),
  set_prior("student_t(3, 0, 0.1)", class = "sd", dpar = "muother")
)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = chains, 
           refresh = refresh,
           prior = priors,
           iter = ifelse(is_ci(), 400, 2400),
           file = here("results/models/Class_year__year_detool_year.rds"),
           file_refit = "on_change"
)

yrep_fit_posterior <- posterior_predict(mod, ndraws = 10)
sum(is.na(yrep_fit_posterior))
length(yrep_fit_posterior)

conditions <- make_conditions(data, vars = "de_tool")
rownames(conditions) <- conditions$de_tool
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions,
                              categorical = TRUE,
                              re_formula = NULL),
          plot = FALSE)
p3a <- p$`year:cats__` + 
  scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(range(pvalues_sample$year)[1], range(pvalues_sample$year)[2], 2)) +
  facet_wrap(~de_tool, nrow = 1) +
  labs(y = "Proportion",
       x = "Year") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
p3a$layers[[1]]$aes_params$alpha <- 0.2

#'
#'
#+
f <- Class ~ de_tool
family <- categorical()
data <- pvalues_sample %>% filter(year >= 2018) # keep only values from 2018 and up!
# get_prior(f, data, family)
priors <- c(
  set_prior("normal(0, 1)", class = "b", dpar = "muuniform")
)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           prior = priors,
           chains = chains, 
           cores = chains, 
           refresh = refresh,
           control = list(adapt_delta = 0.99),
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/Class_detool_2018-19.rds"),
           file_refit = "on_change"
)
p <- plot(conditional_effects(mod, 
                              categorical = TRUE, 
                              effects = "de_tool"), 
          plot = FALSE)
p3b <- p$`de_tool:cats__` + 
  labs(y = "Proportion") +
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.position = "bottom")
p3b$layers[[1]]$aes_params$size <- 1
p3 <- p3a / p3b + plot_annotation(tag_levels = "A") +  plot_layout(guides = 'auto')
ggsave(here("figures/figure_3.pdf"), plot = p3, width = 18, height = 12, units = "cm", dpi = 300)
ggsave(here("figures/figure_3.eps"), plot = p3, width = 18, height = 12, units = "cm", dpi = 300)
ggsave(here("figures/figure_3.png"), plot = p3, width = 18, height = 12, units = "cm", dpi = 300)

#'
#'
#+
f <- pi0 ~ de_tool
family <- student()
data <- pvalues_sample
priors <- set_prior("normal(0, 2)", class="b")
mod <- brm(formula = f, 
           data = data, 
           family = family,
           prior = priors,
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0_detool_sample.rds"),
           file_refit = "on_change"
)
p <- plot(conditional_effects(mod, 
                              effects = "de_tool",
                              re_formula = NULL),
          plot = FALSE)
p4b <- p$de_tool +
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank())
p4b$layers[[1]]$aes_params$size <- 1

#+
p4a <- pvalues_sample %>% 
  filter(Class %in% c("anti-conservative", "uniform")) %>% 
  ggplot() + 
  geom_histogram(aes(pi0), color = "white", binwidth = 0.1) +
  labs(x = expression(pi * 0), y = "Count")

cancer <- read_csv(here("results/data/cancer.csv")) %>% 
  mutate(cancer = 1)
data <- pvalues_sample %>% 
  left_join(cancer) %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer)) %>% 
  select(Accession, id, pi0, cancer) %>% 
  drop_na()

p4c <- data %>% 
  ggplot(aes(pi0, factor(cancer))) +
  stat_histinterval(breaks = 10, point_interval = mean_qi) +
  labs(x = expression(pi * 0)) +
  scale_y_discrete(labels = c("non-cancer", "cancer")) +
  theme(axis.title.y = element_blank())

data %>% 
  group_by(cancer) %>% 
  summarise_at("pi0", list(mean = mean, sd = sd))

#+ Fig4, fig.cap=""
p4 <- p4a + p4b + p4c + plot_annotation(tag_levels = "A")
ggsave(here("figures/figure_4.pdf"), plot = p4, width = 24, height = 8, units = "cm", dpi = 300)
ggsave(here("figures/figure_4.eps"), plot = p4, width = 24, height = 8, units = "cm", dpi = 300)
ggsave(here("figures/figure_4.png"), plot = p4, width = 24, height = 8, units = "cm", dpi = 300)

#' 
#' ## Figure 5 
#+
parsed_suppfiles2 <- parsed_suppfiles_raw %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

pvalues_sample2 <- pvalues %>% 
  select(Accession, id, Set, raw = Class, de_tool) %>% 
  inner_join(
    pvalues_filtered %>% 
      select(Accession, id, Set, filtered = Class, de_tool)
  ) %>% 
  inner_join(
    pvalues_sample %>% 
      select(Accession, id, Set)
  )

make_sankey <- function(data) {
  links <-  data %>% 
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
  
  sankeyNetwork(Links = links, 
                Nodes = nodes, 
                Source = "source",
                Target = "target", 
                Value = "value", 
                NodeID = "name", 
                fontSize = 24, 
                nodeWidth = 30,
                nodePadding = 10,
                margin = list(top = 0, right = 0, bottom = 0, left = 0),
                fontFamily = "Helvetica")
}

#' figure_5A.png
#+
save_sankey_as_webshot <- function(p, path) {
  html <- tempfile(fileext = ".html")
  saveNetwork(p, html)
  webshot(html, path)
}

save_sankey_as_webshot(pvalues_sample2 %>% make_sankey(), here("figures/figure_5A.png"))

get_props <- function(data) {
  raw <- data %>% 
    rename(Class = raw) %>% 
    count(Class, name = "raw")
  filt <- data %>% 
    rename(Class = filtered) %>% 
    count(Class, name = "filtered")
  raw %>% 
    left_join(filt) %>% 
    mutate(prop_raw = raw / sum(raw),
           prop_filt = filtered / sum(filtered),
           FC = filtered / raw)
}

if (!exists("snakemake")) {
  get_props(pvalues_sample2)
}

by_detool <- pvalues_sample2 %>% 
  group_by(de_tool) %>% 
  nest() %>% 
  mutate(sankey = map(data, make_sankey),
         props = map(data, get_props))

by_detool$path <- c("figures/figure_5B.png", "figures/figure_5C.png", "figures/figure_5E.png", "figures/figure_5F.png", "figures/figure_5D.png")
by_detool %>% 
  mutate(res = map2(sankey, path, ~save_sankey_as_webshot(.x, here(.y))))

imgs <- glue::glue("figures/figure_5{LETTERS[1:6]}.png")
png_to_ggplot <- function(path) {
  img <- image_read(here(path))
  ggplot() + 
    annotation_custom(rasterGrob(img)) +
    theme(axis.line = element_blank())
}

plots <- imgs %>% 
  map(png_to_ggplot)

titles <- as.list(c("Total", by_detool %>% arrange(path) %>% pull(de_tool)))
names(titles) <- LETTERS[1:6]

plots2 <- map2(plots, titles, ~ .x + labs(title = .y) + theme(plot.title = element_text(hjust = 0.5, vjust=-6)))
patchwork <- wrap_plots(plots2) + 
  plot_annotation(tag_levels = "A")
ggsave(here("figures/figure_5.pdf"), plot = patchwork, width = 18, height = 12, units = "cm", dpi = 300)
ggsave(here("figures/figure_5.eps"), plot = patchwork, width = 18, height = 12, units = "cm", dpi = 300)
ggsave(here("figures/figure_5.png"), plot = patchwork, width = 18, height = 12, units = "cm", dpi = 300)

rescue_efficiency <- by_detool %>% 
  select(de_tool, props) %>% 
  unnest(props)

#' Figure 6
#' Importing publication data.
#+
publications <- read_csv(
  here("results/data/publications.csv"),
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

#' Importing citation data.
#+
citations <- read_csv(
  here("results/data/scopus_citedbycount.csv"),
  col_types = cols(
    PubMedIds = col_character(),
    citations = col_double()
  ))

#' Import Citescore data.
#+
path <- here("resources/data/CiteScore 2011-2019 new methodology - October 2020.xlsx")
citescore <- readxl::excel_sheets(path) %>% 
  str_subset("^CiteScore")
citescores_raw <- citescore %>% 
  map(~readxl::read_xlsx(path, sheet = .x)) %>% 
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
  filter(year <= 2020)

pubs_citescore <- pubs_all %>% 
  drop_na() %>% 
  left_join(citations)

#' 
#' Importing document summaries.
#+
document_summaries <- read_csv(
  here("results/data/document_summaries.csv"),
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
         year = lubridate::year(PDAT)) %>% 
  group_by(year) %>% 
  summarise(published = sum(published), total = n())

#' Add Accession to PubMedIds
#+
acc_pmid <- document_summaries %>% 
  select(Accession, PubMedIds, taxon) %>% 
  mutate(PubMedIds = map(PubMedIds, str_split, ";"),
         PubMedIds = map(PubMedIds, 1)) %>% 
  unnest(PubMedIds) %>% 
  drop_na()

data <- pvalues_sample %>% 
  select(Accession, Class) %>% 
  left_join(acc_pmid %>% 
              left_join(pubs_citescore) %>% 
              drop_na()) %>% 
  drop_na() %>% 
  mutate(
    anticons = as.numeric(Class == "anti-conservative"),
    year = year - min(year)
  )

f <- anticons ~ CiteScore + year
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
  file = here("results/models/anticons__CiteScore_year.rds"),
  file_refit = "on_change"
)

p <- plot(
  conditional_effects(
    mod, 
    effects = "CiteScore"), 
  line_args = list(color = "black"), plot = FALSE)
pa <- p$CiteScore + 
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Journal CiteScore")

h1 <- hypothesis(mod, "CiteScore > 0")
citescore_es <- data %>% 
  data_grid(CiteScore = c(0.1, 60), year = 0) %>% 
  add_epred_draws(mod) %>% 
  ungroup() %>% 
  select(CiteScore, .draw, .epred) %>% 
  pivot_wider(names_from = CiteScore, values_from = .epred) %>% 
  unnest(cols = c(`0.1`, `60`)) %>% 
  transmute(es = `0.1` - `60`)
citescore_es %>% 
  median_qi()
citescore_es %>% 
  transmute(es = es > 0) %>% 
  mean_qi()
citescore_es %>% 
  transmute(es = es > 0.05) %>% 
  mean_qi()
citescore_es %>% 
  transmute(es = es < -0.05) %>% 
  mean_qi() %>% 
  mutate_at(c("es", ".lower", ".upper"), prettyNum)
citescore_es %>% 
  ggplot() +
  geom_density(aes(es)) +
  geom_vline(xintercept = 0, linetype = "dashed")

f <- anticons ~ log_citations + year
log_citations <- data %>% 
  mutate_at("citations", list(log_citations = ~log10(.x + 0.1))) %>% 
  mutate(year = year - min(year))
mod <- brm(
  formula = f, 
  data = log_citations,
  family = family, 
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = here("results/models/anticons__log_citations_year.rds"),
  file_refit = "on_change"
)
p <- plot(
  conditional_effects(
    mod, 
    effects = "log_citations"), 
  line_args = list(color = "black"), plot = FALSE)
pb <- p$log_citations +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), labels = c(0, 1, 10, 100, 1000, 10000)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Article citations")
o <- hypothesis(mod, "log_citations > 0")

p <- pa + pb + plot_annotation(tag_levels = "A")
ggsave(here("figures/figure_6.pdf"), plot = p, height = 8, width = 12, dpi = 300, units = "cm")
ggsave(here("figures/figure_6.eps"), plot = p, height = 8, width = 12, dpi = 300, units = "cm")
ggsave(here("figures/figure_6.png"), plot = p, height = 8, width = 12, dpi = 300, units = "cm")
