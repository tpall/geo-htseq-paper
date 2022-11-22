library(tidyverse)

parsed_suppfiles <- read_csv("output/parsed_suppfiles.csv")
parsed_suppfiles
raw <- parsed_suppfiles %>% 
  filter(is.na(note), Type == "raw")
raw_acc <- raw %>% 
  mutate(geo_accession = str_extract(id, "GSE\\d+")) %>% 
  count(geo_accession, name = "n_pvalue_sets")
spots <- read_csv("results/data/spots.csv")
d <- spots %>% 
  inner_join(raw_acc) %>% 
  count(geo_accession, name = "n_samples") %>% 
  left_join(raw_acc) %>% 
  mutate(N = n_samples / (n_pvalue_sets + 1)) %>% 
  rename(Accession = geo_accession)

spots_missing <- spots %>% 
  add_count(geo_accession, name = "n_samples") %>% 
  filter(is.na(run_accession)) %>% 
  pull(geo_accession)

document_summaries <- read_csv("output/document_summaries.csv")
document_summaries %>% 
  filter(Accession %in% spots_missing) %>% 
  write_csv("results/document_summaries_spots_rerun.csv")

de_tools_from_summaries <- function(x) {
  x %>% 
  map(str_to_lower) %>% 
  map(str_extract_all, "deseq2?|de(g|x)seq|rockhopper|cuff(diff?|links)|edger|clc ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth") %>% 
  map(unlist) %>% 
  map(unique)
  }

ds1 <- document_summaries %>% 
  select(Accession, taxon, summary) %>% 
  mutate(de_tools = map(summary, de_tools_from_summaries))


ds1 %>% 
  select(Accession, de_tools) %>% 
  unnest(de_tools) %>% 
  unnest(de_tools) %>% 
  group_by(Accession) %>% 
  summarise(de_tools = str_c(de_tools, collapse = ", ")) %>% 
  ungroup() %>% 
  count(de_tools) %>% 
  arrange(desc(n)) %>% 
  view()



de_tools %>% 
  unlist() %>% 
  table()

acc_with_de_tools <- document_summaries %>% 
  filter(map_lgl(de_tools, ~length(.x) > 0)) %>% 
  pull(Accession)

raw %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  filter(Accession %in% acc_with_de_tools)



sample_titles <- document_summaries %>% 
  # head() %>% 
  mutate(titles = str_split(Samples, ",")) %>% 
  unnest(titles) %>% 
  select(Accession, titles, n_samples) %>% 
  filter(str_detect(titles, "Title")) %>% 
  mutate_at("titles", str_remove_all, '"') %>% 
  mutate_at("titles", str_remove_all, ' Title: ') %>% 
  mutate_at("titles", str_remove_all, '[\\}\\]]')

ns <- sample_titles %>% 
  filter(Accession %in% pvalues_sample$Accession) %>% 
  mutate(titles2 = str_remove(titles, "\\d+$"))
ns2 <- ns %>% 
  group_by(Accession, n_samples) %>% 
  mutate(n_unique_titles = n_distinct(titles2),
         N = n_samples / n_unique_titles
         )

ns2 %>% 
  write_csv("results/another_N-table_for_ylo2.csv")

ns2 %>% 
  select(Accession, n_unique_titles, N) %>% 
  distinct() %>% 
  filter(N < 50) %>% 
  ggplot() +
  geom_histogram(aes(N), binwidth = 1)

pvalues_sample <- read_csv(here::here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)

get_n_tests <- function(x) {
  x %>% 
    str_remove_all("[\\[\\]]") %>% 
    str_split(", ", simplify = TRUE) %>% 
    as.numeric() %>% 
    sum(na.rm = TRUE)
}

pvalues_sample %>% 
  mutate(
    n_tests = map_dbl(hist, get_n_tests),
    `prop_pvals<FDR` = ifelse(Class == "anti-conservative", FDR_pval / n_tests, NA_real_)
    ) %>% 
  left_join(d) %>% 
  write_csv("results/another_table_for_ylo.csv")


x <- pvalues_sample$hist
head(x)



get_n_tests(x[[1]])


d %>% 
  ggplot() +
  geom_histogram(aes(n_samples), bins = 100)
d %>% 
  ggplot() +
  geom_jitter(aes(n_samples, n), size = 1/3, width = 1/3) +
  scale_x_log10() +
  scale_y_log10()


d %>% 
  mutate(N = ceiling(n_samples / n)) %>% 
  ggplot() +
  geom_histogram(aes(N), binwidth = 1)

d %>% 
  mutate(N = n_samples / (n + 1)) %>% 
  pull(N) %>% 
  table()


d %>% 
  mutate(N = ceiling(n_samples / (n + 1))) %>% 
  ggplot() +
  geom_histogram(aes(N), binwidth = 1)
