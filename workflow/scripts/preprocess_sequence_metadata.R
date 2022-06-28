if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}

library(readr)
library(dplyr)
library(stringr)
library(here)

spots_raw <- read_csv(here("results/data/spots.csv")) %>% 
  rename_all(str_to_lower)

probs <- problems(spots_raw)

if (nrow(probs) > 0) {
  
  spots_raw <- spots_raw %>% 
    mutate_all(as.character)
  
  #' Fix missing columns
  columns_off <- spots_raw[probs$row, ]
  colnames(columns_off) <- c(setdiff(colnames(columns_off), "sample_name"), "sample_name")
  
  columns_ok <- spots_raw[-probs$row, ]
  
  spots <- bind_rows(columns_ok, columns_off)
} else {
  spots <- spots_raw
}

#' Fix cases with missing platform.
#+
seq_platform <- spots %>% 
  mutate_at(vars(platform, model, library_layout, library_source, library_selection, library_strategy, tax_id), str_to_lower) %>% 
  mutate(
    model = ifelse(is.na(model) & platform == "bgiseq-500", platform, model),
    platform = ifelse(platform == "bgiseq-500" & model == "bgiseq-500", "bgiseq", platform),
    model = ifelse(is.na(model) & platform == "hiseq x ten", platform, model),
    platform = ifelse(platform == "hiseq x ten" & model == "hiseq x ten", "illumina", platform),
    model = ifelse(is.na(model) & str_detect(platform, "illumina"), platform, model),
    platform = ifelse(str_detect(model, "illumina"), "illumina", platform),
    library_layout = ifelse(library_layout %in% c("paired", "single"), library_layout, NA_character_),
    library_layout = ifelse(is.na(library_layout) & library_selection %in% c("paired", "single"), library_selection, library_layout),
    library_selection = ifelse(library_selection %in% c("paired", "single"), library_source, library_selection),
    library_strategy = ifelse(library_strategy == "transcriptomic", "rna-seq", library_strategy),
    tax_id = str_trim(str_remove(tax_id, "\\.0$")),
    tax_id =  case_when(
      tax_id == "mus musculus" ~ "10090",
      tax_id == "acipenser sinensis" ~ "61970",
      tax_id == "aegilops longissima" ~ "4486",
      tax_id == "aegilops tauschii" ~ "37682",
      tax_id == "amycolatopsis mediterranei" ~ "33910",
      tax_id == "amycolatopsis mediterranei u32" ~ "749927",
      tax_id == "aquilegia coerulea" ~ "218851",
      tax_id == "arabidopsis thaliana" ~ "3702",
      tax_id == "aspergillus niger" ~ "5061",
      tax_id == "camellia sinensis" ~ "4442",
      tax_id == "capra hircus" ~ "9925",
      tax_id == "clostridioides difficile r20291" ~ "645463",
      tax_id == "danio rerio" ~ "7955",
      tax_id == "diospyros kaki" ~ "35925",
      tax_id == "drosophila melanogaster" ~ "7227",
      tax_id == "dunaliella salina" ~ "3046",
      tax_id == "escherichia coli" ~ "562",
      tax_id == "ettlia oleoabundans" ~ "1127754",
      tax_id == "gluconobacter oxydans" ~ "442",
      tax_id == "glycine max" ~ "3847",
      tax_id == "haemophilus influenzae 86-028np" ~ "281310",
      tax_id == "halicephalobus mephisto" ~ "2559892",
      tax_id == "homo sapiens" ~ "9606",
      tax_id == "human betaherpesvirus 6a" ~ "32603",
      tax_id == "hymenophyllum caudiculatum" ~ "295381",
      tax_id == "hymenophyllum dentatum" ~ "638559",
      tax_id == "leishmania donovani" ~ "5661",
      tax_id == "macaca fascicularis" ~ "9541",
      tax_id == "macaca mulatta" ~ "9544",
      tax_id == "malus domestica" ~ "3750",
      tax_id == "manis javanica" ~ "9974",
      tax_id == "montastraea cavernosa" ~ "63558",
      tax_id == "mus musculus" ~ "10090",
      tax_id == "oncorhynchus gorbuscha" ~ "8017",
      tax_id == "oryza sativa" ~ "4530",
      tax_id == "oryza sativa japonica group" ~ "39947",
      tax_id == "paeonia lactiflora" ~ "35924",
      tax_id == "papaver somniferum" ~ "3469",
      tax_id == "papilio machaon" ~ "76193",
      tax_id == "papilio xuthus" ~ "66420",
      tax_id == "plasmodium falciparum" ~ "5833",
      tax_id == "plasmodium falciparum 3d7" ~ "36329",
      tax_id == "pleurotus tuoliensis" ~ "879823",
      tax_id == "pseudomonas aeruginosa" ~ "287",
      tax_id == "pseudomonas cichorii" ~ "36746",
      tax_id == "saccharomyces cerevisiae" ~ "4932",
      tax_id == "saccharomyces mikatae" ~ "114525",
      tax_id == "saccharomyces paradoxus" ~ "27291",
      tax_id == "sclerotinia sclerotiorum" ~ "5180",
      tax_id == "solanum lycopersicum" ~ "4081",
      tax_id == "solanum microdontum" ~ "73574",
      tax_id == "solanum microdontum subsp. gigantophyllum" ~ "710637",
      tax_id == "sporothrix schenckii 1099-18" ~ "1397361",
      tax_id == "staphylococcus aureus" ~ "1280",
      tax_id == "tigriopus californicus" ~ "6832",
      tax_id == "triticum"  ~ "4564",
      tax_id == "triticum urartu" ~ "4572",
      tax_id == "vitis vinifera" ~ "29760",
      tax_id == "zea mays" ~ "4577",
      tax_id == "zebrafish metagenome" ~ "7955",
      TRUE ~ tax_id
    )
    )

sequencing_metadata <- seq_platform %>% 
  mutate_at("spots", as.numeric) %>% 
  group_by(geo_accession, library_strategy, library_source, library_selection, library_layout, platform, model, tax_id) %>% 
  summarise_at("spots", list(reads = mean)) %>% 
  ungroup()

sequencing_metadata %>% 
  write_csv(here("results/sequencing_metadata.csv"))

sequencing_metadata_unique_platform <- sequencing_metadata %>% 
  group_by(Accession = geo_accession) %>% 
  add_count() %>% 
  ungroup() %>% 
  filter(n == 1) %>% 
  select(Accession, everything()) %>% 
  select(-geo_accession, -n)

write_csv(sequencing_metadata_unique_platform, here("results/sequencing_metadata_unique_platform.csv"))
