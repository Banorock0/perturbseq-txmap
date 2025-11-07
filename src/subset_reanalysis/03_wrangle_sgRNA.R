#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(janitor)

setwd("~/qb25-answers/perturbseq-txmap/src/subset_reanalysis/")
sgRNA_raw_csv <- read.xlsx(
  "crispri_library.xlsx",
  startRow = 1
) %>% 
  clean_names()


feature_ref_df_a <- sgRNA_raw_csv %>% 
  select(unique_sg_rna_pair_id)


