
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
library(forcats)
library(here)
font_size <- 8
font_family <- "Arial"
old <- theme_set(theme_cowplot(font_size = font_size, font_family = font_family))

if (!is_phantomjs_installed()) {
  install_phantomjs()
} 

#+ params
is_ci <- function() {
  "CI" %in% Sys.getenv()
}
chains <- 3 #ifelse(is_ci(), 1, 4)
cores <- chains
refresh = 200
rstan_options(auto_write = TRUE, javascript = FALSE)
if (!dir.exists("results/models")) {
  dir.create("results/models", recursive = TRUE)
}

#+ data
parse_detool <- . %>% 
  mutate(
    de_tool = case_when(
      is.na(de_tool) ~ "unknown",
      de_tool %in% c("biojupie", "exdega", "nascent rna seq") ~ "unknown",
      de_tool == "cufflinks-edger" ~ "edger",
      de_tool == "cufflinks-deseq" ~ "deseq",
      de_tool == "cufflinks-deseq2" ~ "deseq2",
      de_tool %in% c("cufflinks-limma", "edger-limma") ~ "limma",
      de_tool == "clc genomics-edger" ~ "clc genomics",
      TRUE ~ de_tool
    )
  )

pvalues <- read_csv(here("results/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform) %>% 
  parse_detool()
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform) %>% 
  parse_detool()
pvalues_filtered <- read_csv(here("results/pvalues_filtered.csv")) %>% 
  rename(de_tool = analysis_platform) %>% 
  parse_detool()
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))

#' 
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
  filter(str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())

#+ fig1
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

class_counts <- pvalues_sample %>% 
  count(Class, name = "N")

fit <- brm(Class ~ 1, 
           data = pvalues_sample, 
           family = categorical(), 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/Class_1.rds"),
           file_refit = "never"
)

pe <- posterior_epred(fit, ndraws = 100)
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

fig1a <- tableGrob(
  arrange(plot_data, Class), 
  rows = NULL, 
  theme = ttheme_minimal(base_size = font_size))
fig1a <- gtable_add_grob(fig1a,
                         grobs = segmentsGrob( # line across the bottom
                           x0 = unit(0,"npc"),
                           y0 = unit(0,"npc"),
                           x1 = unit(1,"npc"),
                           y1 = unit(0,"npc"),
                           gp = gpar(lwd = 2.0)),
                         t = 1, b = 1, l = 1, r = 3)

fig1b <- hist_data_plots %>% 
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
  theme_void() +
  theme(strip.text = element_text(size = font_size))

fig1 <- (wrap_elements(fig1b) | (plot_spacer() / wrap_elements(fig1a) / plot_spacer())) + 
  plot_layout(widths = c(2, 3)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(hjust = 1, vjust = 0))
# ggsave(here("figures/Fig1.pdf"), plot = fig1, width = 20, height = 12, units = "cm", dpi = 300)
# ggsave(here("figures/Fig1.eps"), plot = fig1, width = 20, height = 12, units = "cm", dpi = 300)
tiff(here("figures/Fig1.tiff"), width = 1600, height = 1100, units = "px", res = 300)
plot(fig1)
dev.off()

#### Fig. 2 Class~detool
#'
#+ fig2
# We will rescale years, so that year = 0 is the year of first submission
mod <- read_rds(here("results/models/Class_year__year_detool_year.rds"))
draws <- mod$data %>% 
  select(year, de_tool) %>% 
  filter(de_tool %in% c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  distinct() %>% 
  add_epred_draws(mod)
p2a <- draws %>% 
  ggplot(aes(year + min(pvalues_sample$year), .epred)) +
  stat_summary(fun.data = mean_hdci, aes(fill = .category), geom = "ribbon",  alpha = 0.2, fun.args = list(.width = 0.95)) +
  stat_summary(fun = mean, aes(color = .category), geom = "line") +
  facet_wrap(~str_replace(de_tool, "-", "-\n"), nrow = 1, scales = "free_x") +
  labs(y = "Proportion",
       x = "Year") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_colour_manual(values = c("#00BF7D", "#A3A500", "#F8766D", "#00B0F6", "#E76BF3")) +
  scale_fill_manual(values = c("#00BF7D", "#A3A500", "#F8766D", "#00B0F6", "#E76BF3"))

#'
#'
#+
mod <- read_rds(here("results/models/n__trials(total_in_de_tool)__Class_de_tool_Class:de_tool_2018up.rds"))
draws <- mod$data %>% 
  data_grid(Class, de_tool = c("cuffdiff", "deseq", "deseq2", "edger", "limma"), total_in_de_tool = 1000) %>% 
  add_linpred_draws(mod) %>% 
  mutate_at(".linpred", inv_logit_scaled)
p2b <- draws %>% 
  ungroup() %>% 
  filter(Class != "uniform") %>% 
  ggplot(aes(de_tool, .linpred)) +
  stat_pointinterval(point_size = 1) +
  facet_wrap(~Class, nrow = 1) +
  labs(y = "Proportion") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )

p2 <- p2a / p2b + plot_annotation(tag_levels = "A") +  plot_layout(guides = 'auto')
# ggsave(here("figures/Fig2.pdf"), plot = p2, width = 18, height = 12, units = "cm", dpi = 300)
# ggsave(here("figures/Fig2.eps"), plot = p2, width = 18, height = 12, units = "cm", dpi = 300)
tiff(here("figures/Fig2.tiff"), width = 2244, height = 1562, units = "px", res = 300)
plot(p2)
dev.off()

### Fig. 3 pi0
#'
#+ fig3
mod <- read_rds(here("results/models/pi0_detool_sample.rds"))
draws <- mod$data %>% 
  data_grid(de_tool = c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  add_epred_draws(mod)
p3b <- draws %>% 
  ggplot(aes(de_tool, .epred)) +
  stat_pointinterval(point_size = 1) +
  labs(y = expression(pi[0])) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

#+
p3a <- pvalues_sample %>% 
  filter(Class %in% c("anti-conservative", "uniform")) %>% 
  mutate(bins = cut_width(pi0, 0.1, boundary = 0)) %>% 
  count(bins) %>% 
  mutate(p = n / sum(n)) %>% 
  ggplot() + 
  geom_col(aes(bins, p)) +
  labs(x = expression(pi * 0~bins), y = "Proportion of\nP value sets") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

cancer <- read_csv(here("results/data/cancer.csv")) %>% 
  mutate(cancer = 1)
data <- pvalues_sample %>% 
  left_join(cancer) %>% 
  mutate(cancer = ifelse(is.na(cancer), 0, cancer)) %>% 
  select(Accession, id, pi0, cancer) %>% 
  drop_na()

p3c <- data %>% 
  ggplot(aes(pi0, factor(cancer))) +
  stat_histinterval(breaks = 10, point_interval = mean_qi) +
  labs(x = expression(pi * 0)) +
  scale_y_discrete(labels = c("non-\ncancer", "cancer")) +
  theme(axis.title.y = element_blank())

tf <- read_csv(here("results/data/transcription_factor.csv")) %>% 
  mutate(traf = 1)
data <- pvalues_sample %>% 
  left_join(tf) %>% 
  mutate(traf = ifelse(is.na(traf), 0, traf)) %>% 
  select(Accession, id, pi0, traf) %>% 
  drop_na()
p3d <- data %>% 
  ggplot(aes(pi0, factor(traf))) +
  stat_histinterval(breaks = 10, point_interval = mean_qi) +
  labs(x = expression(pi * 0)) +
  scale_y_discrete(labels = c("non-\ntranscription\nfactor", "transcription\nfactor")) +
  theme(axis.title.y = element_blank())

#+
p3 <- (p3a + p3b) / (p3c + p3d) + plot_annotation(tag_levels = "A")
# ggsave(here("figures/Fig3.pdf"), plot = p3, width = 12, height = 10, units = "cm", dpi = 300)
# ggsave(here("figures/Fig3.eps"), plot = p3, width = 12, height = 10, units = "cm", dpi = 300)
tiff(here("figures/Fig3.tiff"), width = 1654, height = 1151, units = "px", res = 300)
plot(p3)
dev.off()

#### Fig. 4 RNA-seq power simulation
n_data <- read_csv(here("results/n_data.csv"))
nplot <- n_data %>% 
  select(Accession, N) %>% 
  distinct() %>% 
  filter(!is.na(N)) %>% 
  mutate(Nbins = fct_lump(factor(N), 10, other_level = ">10")) %>% 
  count(Nbins) %>% 
  ggplot() +
  geom_col(aes(Nbins, n)) +
  labs(x = "N", y = "Count")

simres_df_parsed  <- read_csv(here("results/simres_df_parsed.csv"))
powerplot <- simres_df_parsed %>% 
  ggplot() +
  geom_line(aes(ss1, marginal_power, group = 1 - as.numeric(pde), color = 1 - as.numeric(pde)), size = 0.51) +
  labs(x = "N", y = "Power") +
  scale_color_continuous(expression(True~pi*0)) +
  facet_wrap(~str_to_sentence(set))
power_fig <- nplot + powerplot + plot_annotation(tag_levels = "A")
# ggsave(here("figures/Fig4.pdf"), plot = power_fig, height = 6, width = 12, dpi = 300, units = "cm")
# ggsave(here("figures/Fig4.eps"), plot = power_fig, height = 6, width = 12, dpi = 300, units = "cm")
tiff(here("figures/Fig4.tiff"), width = 1654, height = 781, units = "px", res = 300)
plot(power_fig)
dev.off()

#### Fig. 5 N versus Class
ndata_ac <- n_data %>% 
  filter(Class != "uniform", !is.na(N)) %>%
  mutate(
    Nb = case_when(
      N <= 6 ~ as.character(N), 
      N > 6 & N <= 10 ~ "7-10",
      TRUE ~ ">10"
    ),
    Nb = factor(Nb, levels = c(as.character(1:6), "7-10", ">10"))
  ) %>% 
  select(Accession, N, Nb, anticons) %>% 
  distinct()

nac_mod <- brm(
  anticons ~ Nb,
  data = ndata_ac,
  family = bernoulli(),
  prior = prior("normal(0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/anticons__N.rds"),
  file_refit = "on_change"
)
# pp_check(nac_mod)
# summary(nac_mod)
nac_mod_draws <- nac_mod$data %>% 
  data_grid(Nb) %>% 
  add_epred_draws(nac_mod)
pnacmod_nb <- nac_mod_draws %>% 
  ggplot() +
  stat_pointinterval(aes(Nb, .epred), point_size = 1) +
  labs(x = "N", y = "Proportion of anti-conservative\np value histograms")

nclassb_data <- n_data %>% 
  filter(Class != "uniform", !is.na(N)) %>% 
  mutate(
    Nb = case_when(
      N <= 6 ~ as.character(N), 
      N > 6 & N <= 10 ~ "7-10",
      TRUE ~ ">10"
    ),
    Nb = factor(Nb, levels = c(as.character(1:6), "7-10", ">10"))
  ) %>% 
  select(Accession, Class, Nb) %>% 
  distinct() %>% 
  count(Class, Nb) %>% 
  group_by(Class) %>%
  mutate(nn = sum(n), p = n / nn) 

nclassb_mod <- brm(
  n | trials(nn) ~ Class + Nb + Class:Nb,
  data = nclassb_data,
  family = binomial(),
  prior = prior("normal(0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/n | trials(nn)__Class + Nb + Class:Nb.rds"),
  file_refit = "on_change"
)
# pp_check(nclassb_mod)
nclassb_draws <- nclassb_data %>% 
  data_grid(Class = c("anti-conservative", "bimodal", "other", "conservative"), Nb = c('1', '2', '3', '4', '5', '6', '7-10', '>10'), nn = 600) %>% 
  add_linpred_draws(nclassb_mod) %>% 
  mutate_at(".linpred", inv_logit_scaled)
nclassb <- nclassb_draws %>% 
  ggplot() +
  stat_pointinterval(aes(Nb, .linpred, color = Class), position = position_dodge(0.5), point_size = 1) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7-10", ">10")) +
  labs(x = "N", y = "Proportion") +
  guides(color = guide_legend(nrow = 2)) +
  scale_colour_manual(values = c("#00BF7D", "#A3A500", "#F8766D", "#00B0F6", "#E76BF3")) +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank()
    )


### Fig. 5
nac_fig <- pnacmod_nb + nclassb + 
  plot_layout(widths = c(2, 3)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(hjust = 1.1))
# ggsave(here("figures/Fig5.pdf"), plot = nac_fig, height = 9, width = 22.5, dpi = 300, units = "cm")
# ggsave(here("figures/Fig5.eps"), plot = nac_fig, height = 9, width = 22.5, dpi = 300, units = "cm")
tiff(here("figures/Fig5.tiff"), height = 720, width = 1654, units = "px", res = 300)
plot(nac_fig)
dev.off()

### Fig. 6
parsed_suppfiles_rerun <- read_csv(here("results/data/parsed_suppfiles_rerun.csv")) %>% 
  distinct() %>% 
  filter(Type == "raw") %>% 
  select(-Type)

get_n_tests <- function(x) {
  x %>% 
    str_remove_all("[\\[\\]]") %>% 
    str_split(", ", simplify = TRUE) %>% 
    as.numeric() %>% 
    sum(na.rm = TRUE)
}

pvalues_rerun <- pvalues %>% 
  inner_join(parsed_suppfiles_rerun, by = c("id", "Set"), suffix = c(".orig", ".rerun")) %>% 
  rename_all(str_remove, ".orig") %>% 
  rename(prop_alpha_pval = prop_FDR_pval, alpha_pval = FDR_pval) %>% 
  mutate(
    n_tests = map_dbl(hist.rerun, get_n_tests),
    prop_FDR_pval = ifelse(is.na(pi0), NA_real_, FDR_pval.rerun / n_tests)
  )

pvalues_rerun %>% 
  write_csv(here("results/pvalues_rerun.csv"))


p6a <- simres_df_parsed %>% 
  ggplot() +
  geom_line(aes(1 - pde, pi0, color = ss1, group = ss1)) +
  labs(x = expression(True~pi[0]), y = expression(hat(pi[0]))) +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_continuous("N") +
  facet_wrap(~ str_to_sentence(set))

n_data2 <- pvalues %>% 
  filter(!is.na(pi0)) %>%
  inner_join(n_data %>% select(Accession, id, N)) %>% 
  mutate(
    Nb = case_when(
      N <= 6 ~ as.character(N), 
      N > 6 & N <= 10 ~ "7-10",
      TRUE ~ ">10"
    ),
    Nb = factor(Nb, levels = c(as.character(1:6), "7-10", ">10"))
  )

npi0mod <- brm(
  pi0 ~ Nb,
  data = n_data2,
  family = Beta(),
  prior = prior("normal(0, 1)", class = "b"),
  chains = chains, 
  cores = cores, 
  refresh = refresh,
  iter = ifelse(is_ci(), 400, 2000),
  file = here("results/models/pi0 ~ N.rds"),
  file_refit = "on_change"
)

npi0mod_draw <- n_data2 %>% 
  data_grid(Nb = c('1', '2', '3', '4', '5', '6', '7-10', '>10')) %>% 
  add_epred_draws(npi0mod)

p6b <- npi0mod_draw %>% 
  mutate(Nb = factor(Nb, levels = c('1', '2', '3', '4', '5', '6', '7-10', '>10'))) %>% 
  ggplot() +
  stat_pointinterval(aes(Nb, .epred), point_size = 1) +
  labs(x = "N", y  = expression(pi[0]))

p6 <- (p6a + p6b + plot_layout(widths = c(8, 4), nrow = 1)) + plot_annotation(tag_levels = "A")

# ggsave(here("figures/Fig6.pdf"), plot = p6, height = 9, width = 18, dpi = 300, units = "cm")
# ggsave(here("figures/Fig6.eps"), plot = p6, height = 9, width = 18, dpi = 300, units = "cm")
tiff(here("figures/Fig6.tiff"), height = 652, width = 1654, units = "px", res = 300)
plot(p6)
dev.off()

### Fig. 7 Sankey
pvalues_sample2 <- pvalues %>% 
  select(Accession, id, Set, raw = Class, de_tool) %>% 
  inner_join(
    pvalues_filtered %>% 
      select(Accession, id, Set, filtered = Class, de_tool)
  ) %>% 
  inner_join(
    pvalues_sample %>% 
      select(Accession, id, Set)
  ) %>% 
  mutate(
    de_tool = as.character(de_tool)
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
                fontSize = 46, 
                nodeWidth = 30,
                nodePadding = 10,
                margin = list(top = 0, right = 0, bottom = 0, left = 0),
                fontFamily = "Arial")
}

#' 
#+ 
save_sankey_as_webshot <- function(p, path) {
  html <- tempfile(fileext = ".html")
  saveNetwork(p, html)
  webshot(html, path)
}

save_sankey_as_webshot(pvalues_sample2 %>% make_sankey(), here("figures/Fig7A.png"))

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
  filter(de_tool %in% c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  group_by(de_tool) %>% 
  nest() %>% 
  mutate(sankey = map(data, make_sankey),
         props = map(data, get_props))


by_detool$path <- str_c("figures/Fig7", LETTERS[2:(nrow(by_detool) + 1)], ".png")
by_detool %>% 
  mutate(res = map2(sankey, path, ~save_sankey_as_webshot(.x, here(.y))))

png_to_ggplot <- function(path) {
  img <- image_read(here(path))
  ggplot() + 
    annotation_custom(rasterGrob(img)) +
    theme(axis.line = element_blank())
}

plots <- c("figures/Fig7A.png", by_detool$path) %>% 
  map(png_to_ggplot)

titles <- as.list(c("Total", by_detool %>% arrange(path) %>% pull(de_tool)))
names(titles) <- LETTERS[2:(nrow(by_detool) + 1)]

plots2 <- map2(plots, titles, ~ .x + labs(title = .y) + theme(plot.title = element_text(hjust = 0.5, vjust=-2)))

######## adding effect size plots #######

stopifnot("Missing some model objects, please run scripts/supporting_information.R first!" =
            all(file.exists(here(c("results/models/anticons_detool.rds", 
                                   "results/models/anticons_detool_filtered.rds",
                                   "results/models/pi0_detool_sample.rds",
                                   "results/models/pi0_detool_sample_filtered.rds")
            ))))

anticons_det <- read_rds(here("results/models/anticons_detool.rds"))
anticons_det_filt <- read_rds(here("results/models/anticons_detool_filtered.rds"))
pi0_det <- read_rds(here("results/models/pi0_detool_sample.rds"))
pi0_det_filt <- read_rds(here("results/models/pi0_detool_sample_filtered.rds"))

de_tools <- pvalues_sample %>% 
  select(de_tool) %>% 
  distinct()

draws_anticons_det <- de_tools %>% 
  add_epred_draws(anticons_det, value = "raw") %>% 
  ungroup() %>% 
  select(de_tool, .draw, raw)
draws_anticons_det_filt <- anticons_det_filt$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(anticons_det_filt, value = "filtered") %>% 
  ungroup() %>% 
  select(de_tool, .draw, filtered)
draws_pi0_det <- pi0_det$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(pi0_det, value = "raw") %>% 
  ungroup() %>% 
  select(de_tool, .draw, raw)
draws_pi0_det_filt <- pi0_det_filt$data %>% 
  select(de_tool) %>% 
  distinct() %>% 
  add_epred_draws(pi0_det_filt, value = "filtered") %>% 
  ungroup() %>% 
  select(de_tool, .draw, filtered)

draws_ac_merged <- draws_anticons_det %>% 
  inner_join(draws_anticons_det_filt)

pd <- position_dodge(0.6)
pg <- draws_ac_merged %>% 
  filter(de_tool %in% c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  pivot_longer(cols = c("raw", "filtered")) %>% 
  mutate(name = factor(name, levels = c("raw", "filtered"))) %>% 
  ggplot(aes(value, de_tool)) +
  stat_pointinterval(aes(color = name), point_size = 1, position = pd) +
  labs(x = "Prop. anti-cons.") +
  scale_y_discrete(limits = rev) +
  scale_color_discrete("P value\nset") +
  theme(axis.title.y = element_blank())

ph <- draws_ac_merged %>% 
  filter(de_tool %in% c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  mutate(es = filtered - raw) %>% 
  ggplot(aes(es, de_tool)) +
  stat_pointinterval(point_size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1/3) +
  labs(x = "Effect size") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-0.1, 0.45)) +
  theme(axis.title.y = element_blank())

draws_pi0_merged <- draws_pi0_det %>% 
  inner_join(draws_pi0_det_filt)

pi <- draws_pi0_merged %>% 
  filter(de_tool %in% c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  pivot_longer(cols = c("raw", "filtered")) %>% 
  mutate(name = factor(name, levels = c("raw", "filtered"))) %>% 
  ggplot(aes(value, de_tool)) +
  stat_pointinterval(aes(color = name), point_size = 1, position = pd) +
  labs(x = expression(pi[0])) +
  scale_y_discrete(limits = rev) +
  scale_color_discrete("P value\nset") +
  theme(axis.title.y = element_blank())

pj <- draws_pi0_merged %>% 
  filter(de_tool %in% c("cuffdiff", "deseq", "deseq2", "edger", "limma")) %>% 
  mutate(es = filtered - raw) %>% 
  ggplot(aes(es, de_tool)) +
  stat_pointinterval(point_size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1/3) +
  labs(x = "Effect size") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-0.1, 0.45)) +
  theme(axis.title.y = element_blank())

######## end of effect size plots ######

patchwork <- wrap_plots(plots2, nrow = 2) / ((plot_spacer() + (pg + ph + pi + pj + plot_layout(nrow = 2, guides = "collect")) + plot_spacer()) + plot_layout(widths = c(1/12, 5/6, 1/12))) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(hjust = 1))

# ggsave(here("figures/Fig7.pdf"), plot = patchwork, width = 20, height = 24, units = "cm", dpi = 300)
# ggsave(here("figures/Fig7.eps"), plot = patchwork, width = 20, height = 24, units = "cm", dpi = 300)
tiff(here("figures/Fig7.tiff"), width = 2244, height = 2625, units = "px", res = 300)
plot(patchwork)
dev.off()

rescue_efficiency <- by_detool %>% 
  select(de_tool, props) %>% 
  unnest(props)

