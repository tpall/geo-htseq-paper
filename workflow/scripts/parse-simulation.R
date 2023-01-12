library(tidyverse)
library(PROPER)
library(qvalue)
library(limma)
library(DESeq2)
library(edgeR)
library(here)
library(cowplot)
old <- theme_set(theme_cowplot())

rds <- list.files("results", pattern = "(bottomly|gilad)_0\\.\\d_\\d+_rerun3$", full.names = TRUE)
simres_df <- tibble(rds = rds) %>% 
  mutate(
    set = str_replace(rds, "output\\/(bottomly|gilad)_(0\\.\\d)_(\\d+)_rerun3$", "\\1"),
    pDE = as.numeric(str_replace(rds, "output\\/(bottomly|gilad)_(0\\.\\d)_(\\d+)_rerun3$", "\\2")),
    Nreps = str_replace(rds, "output\\/(bottomly|gilad)_(0\\.\\d)_(\\d+)_rerun3$", "\\3"),
    simres = map(rds, read_rds)
    )

simres_pi0 <- function(x) {
  pv <- x$pvalue
  pvfilt <- pv[is.finite(x$sim.opts$lBaselineExpr), , ]
  mean(sapply(1:20, function(y) propTrueNull(pvfilt[,y])))
}

simres_qvalue <- function(x) {
  pv <- x$pvalue
  pvfilt <- pv[is.finite(x$sim.opts$lBaselineExpr), , ]
  qs <- lapply(1:20, function(y) qvalue(pvfilt[, y], fdr.level = 0.05))
  mean(unlist(sapply(qs, function(x) mean(x$significant))))
}

simres_df_parsed <- simres_df %>% 
  mutate(
    powers = map(simres, comparePower, alpha.type = "fdr", alpha.nominal = 0.05, stratify.by = "expr", delta = 0.5),
    powers_sum = map(powers, summaryPower),
    powers_sum = map(powers_sum, as_tibble),
    pi0 = map_dbl(simres, simres_pi0),
    pfdr = map_dbl(simres, simres_qvalue)
  ) %>% 
  dplyr::select(set, pDE, powers_sum, pi0, pfdr) %>% 
  unnest(powers_sum) %>% 
  rename_all(str_to_lower) %>% 
  rename_all(str_replace_all, "\\s+", "_")

simres_df_parsed %>% 
  write_csv(here("results/simres_df_parsed.csv"))

simres_df_parsed %>% 
  ggplot() +
  geom_line(aes(ss1, marginal_power, group = 1 - as.numeric(pde), color = 1 - as.numeric(pde))) +
  labs(x = "N", y = "Power") +
  scale_color_continuous(expression(True~pi[0])) +
  facet_wrap(~str_to_sentence(set))

simres_df_parsed %>% 
  ggplot() +
  geom_line(aes(1 - pde, pi0, color = ss1, group = ss1)) +
  labs(x = expression(True~pi[0]), y = expression(hat(pi[0]))) +
  geom_abline(slope = 1, linetype = "dashed") +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_continuous("N") +
  facet_wrap(~ str_to_sentence(set))

simres_df_parsed %>% 
  ggplot() +
  geom_line(aes(1 - pde, 1 - pfdr, color = ss1, group = ss1)) +
  labs(x = expression(True~pi[0]), y = expression(Prop.~nonsignificant)) +
  geom_abline(slope = 1, linetype = "dashed") +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_continuous("N") +
  facet_wrap(~ str_to_sentence(set))
