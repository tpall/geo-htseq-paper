library(tidyverse)
library(PROPER)
library(DESeq2)

run_simu <- function(p.DE, set, Nreps, nsims = 2, path = "output") {
  file <- file.path(path, paste(set, p.DE, Nreps, "rerun3", sep = "_"))
  if (file.exists(file)) {
    simres <- read_rds(file)
  } else {
    seed <- sample(1:1e9, 1)
    sim.opts <- RNAseq.SimOptions.2grp(
      ngenes = 10000, 
      p.DE = p.DE,
      lOD = set, 
      lBaselineExpr = set,
      lfc = function(x) rnorm(x, mean = 0, sd = 1.5),
      sim.seed = seed
    )
    simres <- runSims(
      Nreps = Nreps,
      sim.opts = sim.opts,
      DEmethod = "DESeq2", 
      nsims = nsims
    )
    write_rds(simres, file)
  }
  return(simres)
}

pde <- seq(0.1, 0.9, 0.1)
set <- c("gilad", "bottomly")
Nreps <- c(2, 4, 6, 8, 10)
simres_df <- expand.grid(set = set, pde = pde, Nreps = Nreps) %>% 
  mutate(
    set = as.character(set), 
    simres = pmap(list(pde, set, Nreps), ~run_simu(..1, ..2, ..3, nsims = 50, path = "results"))
  )
