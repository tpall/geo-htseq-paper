library("tidyverse")
library("sparkline")
library("formattable")
library("here")

document_summaries <- read_csv("data/document_summaries.csv")
pvals_wide <- read_csv("output/pvals_wide.csv")
spots <- read_csv("data/spots.csv")
