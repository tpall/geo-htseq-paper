if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}

#+ libs
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(lubridate)
library(here)
library(networkD3)
library(webshot)
library(ggplot2)
library(grid)
library(patchwork)
library(magick)
library(cowplot)
old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

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
f5a <- pvalues_sample %>% 
  make_sankey()
saveNetwork(f5a, here("figures", "figure_5A.html"))
webshot(here("figures", "figure_5A.html"), here("figures", "figure_5A.png"))

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

get_props(pvalues_sample)

de_tools <- parsed_suppfiles %>% 
  filter(!is.na(Conversion), Type != "raw") %>% 
  mutate(analysis_platform = case_when(
    Type == "basemean" ~ "deseq",
    Type == "aveexpr" ~ "limma",
    Type == "logcpm" ~ "edger",
    Type == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
    TRUE ~ "unknown"
  )) %>% 
  select(Accession, id, analysis_platform) %>% 
  distinct()

pvalues_sample_detools <- pvalues_sample %>% 
  left_join(de_tools)

by_detool <- pvalues_sample_detools %>% 
  group_by(analysis_platform) %>% 
  nest() %>% 
  mutate(sankey = map(data, make_sankey),
         props = map(data, get_props))
sankey_plots <- by_detool %>% 
  pull(sankey)

by_detool$fig <- c("figure_5C.png", "figure_5D.png", "figure_5B.png", "figure_5E.png", "figure_5F.png")
by_detool$html <- c("figure_5C.html", "figure_5D.html", "figure_5B.html", "figure_5E.html", "figure_5F.html")
by_detool %>% 
  mutate(html_path = here("figures", html),
         fig_path = here("figures", fig),
         save_html = map2(sankey, html_path, saveNetwork),
         save_png = map2(html_path, fig_path, webshot)
         )


imgs <- glue::glue("figures/figure_5{LETTERS[1:6]}.png")
png_to_ggplot <- function(path) {
  img <- image_read(here(path))
  ggplot() + 
    annotation_custom(rasterGrob(img)) +
    theme(axis.line = element_blank())
}
plots <- imgs %>% 
  map(png_to_ggplot)

titles <- list(A = "Total",
               B = "cuffdiff",
               C = "deseq",
               D = "edger",
               E = "limma",
               F = "unknown")

plots2 <- map2(plots, titles, ~ .x + labs(title = .y) + theme(plot.title = element_text(hjust = 0.5, vjust=-6)))
patchwork <- wrap_plots(plots2) + 
  plot_annotation(tag_levels = "A")
ggsave(here("figures/figure_5.pdf"), plot = patchwork, width = 18, height = 12, units = "cm", dpi = 300)

rescue_efficiency <- by_detool %>% 
  select(analysis_platform, props) %>% 
  unnest(props)

#' How much p values were lost by filtering
#+
parse_hist <- function(x) {
  x <- str_remove_all(x, "[:punct:]")
  x <- str_split(x, pattern = " ")[[1]]
  sum(as.numeric(x))
}

drop_out <- parsed_suppfiles %>% 
  filter(!is.na(Conversion)) %>% 
  select(Accession, id, Type, Set, hist) %>% 
  mutate(n = map_dbl(hist, parse_hist),
         raw_filt = if_else(Type == "raw", "raw", "filtered")) %>% 
  select(-hist)

set.seed(11)
drop_out_sample <- drop_out %>% 
  group_by(Accession, id) %>% 
  mutate(expval = Type[[2]]) %>% 
  select(-Type) %>% 
  pivot_wider(names_from = raw_filt, values_from = n) %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  mutate(prop_out = 1 - filtered/raw)

drop_out_sample <- drop_out_sample %>% 
  mutate(analysis_platform = case_when(
    expval == "basemean" ~ "deseq",
    expval == "aveexpr" ~ "limma",
    expval == "logcpm" ~ "edger",
    expval == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
    TRUE ~ "unknown"
  ))

drop_out_sample %>% 
  ggplot() +
  geom_histogram(aes(prop_out)) +
  facet_wrap(~analysis_platform, scales = "free_y")
