#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(janitor)


# for own use only
setwd("~/perturbseq-txmap/src/subset_reanalysis/")

# load in the sgRNA sequences
sgRNA_raw_csv <- read.xlsx(
  "crispri_library.xlsx",
  startRow = 1
) %>% 
  clean_names()

# get protospacer 1
# protospacer 1 is followed by constant region 3, based of off pJR89 and pJR101
# construct:
# mU6-protospacer_1-CR3
# hU6-protospacer_2-CR1

feature_ref_df_a <- sgRNA_raw_csv %>% 
  select(unique_sg_rna_pair_id, gene, ensembl_gene_id, sg_id_a, targeting_sequence_a) %>% 
  mutate(
    # add CR3 as pattern
    pattern = "GTTTCAGAGCTAAGCACAAG"
  )

sg_a_intermediate <- feature_ref_df_a %>% 
  mutate(
    id = paste0(gene, "-1"), # -1 designates protospacer a
    name = paste0(gene,"-1"), 
    read = "R2", # protospacer sequenced on read 2,
    pattern = paste0("(BC)", pattern), # (BC) refers to protospacer, pattern is constant region 3 (AKA tracrRNA)
    sequence = targeting_sequence_a, # protospacer
    feature_type = "CRISPR Guide Capture",
    target_gene_id = ensembl_gene_id, # ensembl gene ID to ID perturbation
    target_gene_name = gene # common name to ID perturbation
  ) %>% 
  select(id:read, pattern, sequence:target_gene_name)

# get protospacer 2

feature_ref_df_b <- sgRNA_raw_csv %>% 
  select(
    unique_sg_rna_pair_id,
    gene, 
    ensembl_gene_id,
    sg_id_b,
    targeting_sequence_b
  ) %>% 
  mutate(
    pattern = "GTTTAAGAGCTAAGCTGGAA"
  )

sg_b_intermediate <- feature_ref_df_b %>% 
  mutate(
    id = paste0(gene, "-2"), # -1 designates protospacer a
    name = paste0(gene,"-2"), 
    read = "R2", # protospacer sequenced on read 2,
    pattern = paste0("(BC)", pattern), # (BC) refers to protospacer, pattern is constant region 3 (AKA tracrRNA)
    sequence = targeting_sequence_b, # protospacer
    feature_type = "CRISPR Guide Capture",
    target_gene_id = ensembl_gene_id, # ensembl gene ID to ID perturbation
    target_gene_name = gene # common name to ID perturbation
  ) %>% 
  select(id:read, pattern, sequence:target_gene_name)


# generate feature-reference CSV for cell ranger to 
# ID sgRNA libarries

feature_ref_csv <- rbind(sg_a_intermediate, sg_b_intermediate)

write_csv(
  feature_ref_csv,
  paste0(
    format(Sys.time(), "%Y%m%d"),
    "_dualguide_sgRNA_feature_ref.csv"
  )
)
