library(tidyverse)
library(lubridate)
library(brms)
library(here)

theme_set(theme_classic(base_size = 8))

#' Number of unique GEOs
#+
document_summaries <- read_csv("data/document_summaries.csv")
document_summaries %>% 
  filter(PDAT<="2019-12-31") %>% 
  pull(Accession) %>% 
  n_distinct()

document_summaries %>% 
  filter(PDAT<="2019-12-31") %>% 
  group_by(year(PDAT)) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(cumsum = cumsum(n),
         perc = n / cumsum)

#' Supplementary files. We will throw out also all 'filelist.txt' files.
#+
suppfilenames <- read_lines("data/suppfilenames.txt") %>% 
  tibble(suppfilenames = .) %>% 
  filter(
    str_detect(suppfilenames, "suppl"),
    !str_detect(suppfilenames, "suppl/filelist.txt")
  )

acc_year <- document_summaries %>% 
  select(Accession, PDAT) %>% 
  mutate(year = year(PDAT))

#' We right join to keep only subset from our time frame.
conformity <- suppfilenames %>% 
  mutate(Accession = str_to_upper(str_extract(suppfilenames, "GS[Ee]\\d+"))) %>% 
  right_join(acc_year) %>%
  mutate(
    conforms = !(is.na(suppfilenames) | str_detect(str_to_lower(suppfilenames), "readme|_raw.tar$|\\.[bs]am$|\\.bed$|\\.fa(sta)?"))
  )

#' Total number of files conforming
sum(conformity$conforms)

#' Total number of GEOs conforming
conformity %>% 
  filter(conforms) %>% 
  pull(Accession) %>% 
  n_distinct()

#' Number of GEOs conforming per year
conformity_acc <- conformity %>% 
  group_by(Accession, year) %>% 
  summarise(
    conforms = case_when(
      any(conforms) ~ 1,
      TRUE ~ 0
  )) %>% 
  ungroup()
write_csv(conformity_acc, here("output/conformity_acc.csv"))

#' Total number of conforming GEO accessions.
#+
total_conforming <- conformity_acc %>% 
  count(conforms)

#' Conforming GEO submissions per year.
#+
conformity_acc %>% 
  group_by(year) %>% 
  summarise(conforms = sum(conforms),
            n = n(),
            perc = conforms / n)

#' Figure 1
#' Model in scripts/summary_stats_models.R
fit <- brm(conforms ~ year, data = conformity_acc, family = bernoulli(), file = "models/conforms_year.rds")
p <- plot(conditional_effects(fit), plot = FALSE)$year
p + 
  geom_smooth(color = "black") +
  labs(x = "Year", y = "Proportion of submissions conforming\nwith GEO submission guidelines") +
  scale_x_continuous(breaks = seq(2006, 2019, by = 2)) +
  scale_y_continuous(limits = c(0, 1))
ggsave("plots/conforming_per_year.png", height = 7, width = 11, dpi = 300, units = "cm")


#' Number of sets with p-values, 
parsed_suppfiles_raw <- read_csv(here("data/parsed_suppfiles.csv"))
parsed_suppfiles <- parsed_suppfiles_raw %>% 
  filter(!str_detect(id, "_RAW.tar")) %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

#' Parse analysis platform
get_var <- function(x) {
  NAs <- is.na(x)
  if (all(NAs)) {
    return("unknown")
  }
  return(names(x)[!NAs])
}

pi0 <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  select(Accession, id, Set, pi0) %>% 
  na.omit()

hist <- parsed_suppfiles %>% 
  filter(Type == "raw") %>% 
  select(Accession, id, Set, FDR_pval, hist)
  
pvalues <- parsed_suppfiles %>% 
  filter(!is.na(Type)) %>% 
  select(Accession, id, Set, Type, Class) %>%
  pivot_wider(names_from = Type, values_from = Class) %>% 
  group_by(id, Set) %>% 
  mutate(var = get_var(c("logcpm"=logcpm, "rpkm"=rpkm, "fpkm"=fpkm, "basemean"=basemean, "aveexpr"=aveexpr))) %>%
  ungroup() %>% 
  mutate(analysis_platform = case_when(
    var == "basemean" ~ "deseq",
    var == "aveexpr" ~ "limma",
    var == "logcpm" ~ "edger",
    var == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
    TRUE ~ "unknown"
  )) %>% 
  select(Accession, id, Set, Class = raw, analysis_platform) %>% 
  left_join(acc_year) %>% 
  left_join(pi0) %>% 
  left_join(hist) %>% 
  select(-PDAT) %>% 
  mutate(anticons = case_when(
    Class == "anti-conservative" ~ 1,
    TRUE ~ 0
  ))
write_csv(pvalues, "output/pvalues.csv")

#' Number of unique GEO ids imported
geo_import <- parsed_suppfiles %>% 
  filter(is.na(Type) | str_detect(Type, "raw"))
geo_import %>% 
  pull(Accession) %>% 
  n_distinct()

#' Number of unique files downloaded
geo_import %>% 
  select(Accession, id) %>% 
  mutate(file = str_remove(id, "^.+ from ")) %>% 
  pull(file) %>% 
  n_distinct()

#' Number of unique files imported
geo_import %>% 
  pull(id) %>% 
  n_distinct()


#' notes
geo_import %>% 
  count(note) %>% 
  pull(note)

#' Files which we were able to import.
geo_import %>% 
  mutate(imported = case_when(
    str_detect(str_to_lower(note), "error|i\\/o operation|not determine delim|codec can't decode|empty table|missing optional dependency") ~ "fail",
    is.na(note) ~ "yes",
    TRUE ~ "yes"
  )) %>% 
  group_by(Accession, id) %>% 
  summarise(
    imported = case_when(
      any(imported == "yes") ~ 1,
      TRUE ~ 0
    )) %>% 
  ungroup() %>% 
  count(imported)

#' P values
#+

#' Number of distinct GEO-s with p values
n_distinct(pvalues$Accession)

#' Summary stats of p values per GEO
pvalues %>% 
  count(Accession) %>% 
  summarise_at("n", list(mean = mean, median = median, max = max))

#' GEOs with single or 1-3 sets
#+
pvalues %>% 
  count(Accession) %>% 
  count(n) %>% 
  mutate(p = nn / sum(nn),
         ps = cumsum(p))

#' Sample 1 p value set per accession
set.seed(11)
pvalues_sample <- pvalues %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  ungroup()
write_csv(pvalues_sample, here("output/pvalues_sample.csv"))

#' Save table ids for figure
#+
pvalues_sample %>% 
  select(id, Set) %>% 
  write_csv("output/pvalues_acc.csv")

pvalues_sample %>% 
  count(Class) %>% 
  summarise(Class, 
            n,
            p = signif(n / sum(n), 2))
