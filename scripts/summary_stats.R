library(tidyverse)
library(lubridate)

theme_set(theme_classic())

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

#' Supplementary files. We throw out also filelist.txt files.
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

conformity <- suppfilenames %>% 
  mutate(Accession = str_extract(suppfilenames, "GSE\\d+")) %>% 
  full_join(acc_year) %>%
  mutate(
    conforms = !(is.na(suppfilenames) | str_detect(str_to_lower(suppfilenames), "readme|raw.tar|\\.bam$|\\.sam$|\\.bed$|\\.fa(sta)?"))
  )

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
  na.omit() %>% 
  summarise(conforms = sum(conforms),
            n = n(),
            perc = conforms / n)


library(brms)
fit <- brm(conforms ~ year, data = conformity_acc, family = bernoulli())
p <- plot(conditional_effects(fit), plot = FALSE)$year
p + 
  labs(x = "Year", y = "Proportion of submissions conforming\nwith GEO submission guidelines") +
  scale_x_continuous(breaks = seq(2006, 2019, by = 2))

#' Number of sets with p-values
parsed_suppfiles <- read_csv("data/parsed_suppfiles.csv") %>% 
  filter(is.na(Type) | str_detect(Type, "raw")) %>% 
  mutate(Accession = str_extract(id, "GSE\\d+")) %>% 
  select(Accession, everything())
parsed_suppfiles %>% 
  pull(Accession) %>% 
  n_distinct()

parsed_suppfiles %>% 
  pull(Accession) %>% 
  n_distinct()

#' Files which we were able to import.
#+
parsed_suppfiles %>% 
  mutate(imported = case_when(
    str_detect(str_to_lower(note), "error|i\\/o operation|not determine delim|codec can't decode|empty table|missing optional dependency") ~ "fail",
    is.na(note) ~ "yes",
    TRUE ~ "yes"
  )) %>% 
  count(imported)

