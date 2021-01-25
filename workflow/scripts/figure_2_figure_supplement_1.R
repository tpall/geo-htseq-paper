library(tidyverse)

old <- theme_set(theme_cowplot(font_size = 8, font_family = "Helvetica"))

files <- list.files("results/simulation/", pattern = "group1-vs-group2.diffexp.tsv", full.names = TRUE, recursive = TRUE)
de_list <- files %>% 
  map(read_table2, skip = 1, col_names = FALSE)
names(de_list) <- dirname(files) %>% 
  str_split("/") %>% 
  map_chr(4)

de_raw <- de_list %>% 
  bind_rows(.id = "id")

colnames(de_raw) <- c("id", "transcriptid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

de_raw <- de_raw %>% 
  mutate_at("transcriptid", str_remove_all, '\"')

files <- list.files("results/simulation/", pattern = "sim_tx_info.txt", full.names = TRUE, recursive = TRUE)
tx_list <- files %>% 
  map(read_table2, skip = 1, col_names = FALSE)
names(tx_list) <- dirname(files) %>% 
  str_split("/") %>% 
  map_chr(4)
sim <- tx_list %>% 
  bind_rows(.id = "id")
colnames(sim) <- c("id", "transcriptid", "foldchange.1", "foldchange.2", "DEstatus.1", "DEstatus.2")

de <- de_raw %>% 
  separate("id", c("reps", "n_eff"), sep = "_") %>% 
  mutate_at(c("reps", "n_eff"), str_extract, "\\d+") %>% 
  mutate_at(c("reps", "n_eff"), as.numeric) %>% 
  mutate_at("n_eff", ~ 2 * n_eff)


sim <- sim %>% 
  separate("id", c("reps", "n_eff"), sep = "_") %>% 
  mutate_at(c("reps", "n_eff"), str_extract, "\\d+") %>% 
  mutate_at(c("reps", "n_eff"), as.numeric)

qc_threshold <- function(n_pvals, bins, fdr) {
  qbinom(1 - 1 / bins * fdr, n_pvals, 1 / bins)
}

de_plot_data <- de %>% 
  mutate(
    reps = factor(str_c("N=", reps), levels = c("N=3", "N=6", "N=10")),
    n_eff = factor(str_c("Effects=", n_eff), levels = c("Effects=100", "Effects=200", "Effects=400", "Effects=800"))
  )

qc_thr <- qc_threshold(n_pvals = 20000, bins = 40, fdr = 0.05)

de_plot_data %>% 
  ggplot() +
  geom_histogram(aes(pvalue), binwidth = 1/40, center = 1/80) +
  geom_hline(yintercept = qc_thr, color = "red") +
  facet_grid(reps~n_eff) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

de %>% 
  filter(padj <= 0.05) %>% 
  ggplot() +
  geom_histogram(aes(log2FoldChange), bins = 30) +
  facet_grid(reps~n_eff)


de %>% 
  filter(padj <= 0.05) %>% 
  group_by(log2FoldChange < 0, reps, n_eff) %>% 
  summarise_at("log2FoldChange", list(mean = mean, sd = sd))

merged <- de %>% 
  left_join(sim)


merged %>% 
  filter(padj < 0.05) %>% 
  count(reps, n_eff, DEstatus.1, DEstatus.2) %>% 
  mutate(ok = map2_lgl(DEstatus.1, DEstatus.2, any)) %>% 
  group_by(reps, n_eff, ok) %>% 
  summarise_at("n", sum) %>% 
  mutate(p = n / sum(n))

