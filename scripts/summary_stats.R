library(tidyverse)
library(lubridate)
library(brms)

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
    conforms = !(is.na(suppfilenames) | str_detect(str_to_lower(suppfilenames), "readme|_raw.tar$|\\.bam$|\\.sam$|\\.bed$|\\.fa(sta)?"))
  )

#' Total number of files conforming
conformity %>% 
  filter(conforms) %>% 
  nrow()

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
  ))

#' Total number of conforming GEO accessions.
#+
total_conforming <- conformity_acc %>% 
  ungroup() %>% 
  count(conforms)

#' Conforming GEO submissions per year.
#+
conformity_acc %>% 
  group_by(year) %>% 
  summarise(conforms = sum(conforms),
            n = n(),
            perc = conforms / n)


fit <- brm(conforms ~ year, data = conformity_acc, family = bernoulli())
p <- plot(conditional_effects(fit), plot = FALSE)$year
p + 
  geom_smooth(color = "black") +
  labs(x = "Year", y = "Proportion of submissions conforming\nwith GEO submission guidelines") +
  scale_x_continuous(breaks = seq(2006, 2019, by = 2)) +
  scale_y_continuous(limits = c(0, 1))
ggsave("plots/conforming_per_year.png", height = 7, width = 11, dpi = 300, units = "cm")

#' Number of sets with p-values, 
parsed_suppfiles <- read_csv("data/parsed_suppfiles.csv")

raw_sets <- parsed_suppfiles %>% 
  filter(is.na(Type) | str_detect(Type, "raw")) %>% 
  filter(!str_detect(id, "_RAW.tar")) %>% 
  mutate(Accession = str_to_upper(str_extract(id, "GS[Ee]\\d+"))) %>% 
  select(Accession, everything())

#' Number of unique GEO ids imported
raw_sets %>% 
  pull(Accession) %>% 
  n_distinct()

#' Number of unique files imported
raw_sets %>% 
  pull(id) %>% 
  n_distinct()


raw_sets %>% 
  select(Accession, id) %>% 
  filter(str_detect(id, "from"))

raw_sets %>% 
  select(Accession, id) %>% 
  mutate(file = str_remove(id, "^.+ from ")) %>% 
  filter(!str_detect(file, "^GSE"))


#' notes
raw_sets %>% 
  count(note) %>% 
  pull(note)

#' Files which we were able  to import.
raw_sets %>% 
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
pvalues <- raw_sets %>% 
  filter(Type == "raw")

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
  summarise(p = nn / sum(nn),
            ps = cumsum(p))

#' Sample 1 p value set per accession
set.seed(11)
pvalues_acc <- pvalues %>% 
  group_by(Accession) %>% 
  sample_n(1) %>% 
  ungroup()

#' Save table ids for figure
#+
pvalues_acc %>% 
  select(id, Set) %>% 
  write_csv("output/pvalues_acc.csv")

pvalues_acc %>% 
  count(Class) %>% 
  summarise(Class, 
            n,
            p = signif(n / sum(n), 2))

pvalues_acc %>% 
  count(Set) %>% 
  arrange(desc(n))

pvalues_acc_year <- pvalues_acc %>% 
  left_join(acc_year)

fit2 <- brm(Class ~ year, data = pvalues_acc_year, family = categorical())
p <- plot(conditional_effects(fit2, categorical = TRUE), plot = FALSE)$year
p + 
  labs(x = "Year", y = "Proportion of GEO\nsubmissions") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  theme(legend.position = "bottom",
        legend.title = element_blank())
ggsave("plots/class_per_year.png", height = 7, width = 11, dpi = 300, units = "cm")
