#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(janitor)
library(BiocManager)

#BiocManager::install("biomaRt")
#library(biomaRt)

# TODO: Replact select with dplyr::select
# NOTE: Droped guide RNAs with no annontated ensembl ID, (adding them manually crashes cellranger)
# Likely since they were not annotated with the human genome
# used in Repogle 2022

# # get ensembl ids for genes which don't already have them for some reason
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# listFilters(ensembl)
# 
# # test
# getBM(attributes='ensembl_gene_id', 
#       filters = 'hgnc_symbol', 
#       values = "C8orf44", 
#       mart = ensembl)




# for own use only
setwd("~/perturbseq-txmap/src/subset_reanalysis/")

# load in the sgRNA sequences
sgRNA_raw_csv <- read.xlsx(
  "crispri_library.xlsx",
  startRow = 1
) %>% 
  clean_names() %>% 
  filter(either_guide_duplicated == FALSE) %>%  # removes any sgRNA where the guide matches another
  mutate(
    gene = str_extract(
      unique_sg_rna_pair_id, r"((\d+)_([^_]+)_)", group = 2
    ),
    transcript = str_replace_all(
      transcript, ",", "_"
    )
  ) %>% 
  mutate(
    transcript = str_glue(
      "{gene}_{transcript}"
    ),
    gene = str_replace_all(
      gene, ",", "_"
    )
  )

### TEST CODE

# fill in ensembl id if it is missing
# gene_names <- sgRNA_raw_csv$gene
# ensembl_ids <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"), 
#                      filters = 'external_gene_name', 
#                      values = gene_names, 
#                      mart = ensembl) %>% 
#   distinct(ensembl_gene_id, .keep_all= TRUE)
# 
# new_ensembl_ids <- ensembl_ids
# for (i in 1:length(gene_names)) {
#   if (is.na(ensembl_ids[i])) {
#     new_ensembl_ids[i] = getBM(attributes='ensembl_gene_id', 
#                        filters = 'external_gene_name', 
#                        values = gene_names[i], 
#                        mart = ensembl)
#   }
# }

# new_ensembl_ids <- unlist(new_ensembl_ids)
# test <- cbind(gene_names, ensembl_ids, new_ensembl_ids)
# 
# test %>% 
#   as_tibble() %>% 
#   filter(is.na(ensembl_ids))
# 
# # sgRNA_fixed_ensembl <- sgRNA_raw_csv %>% 
# #   mutate(
# #     new_id = new_ensembl_ids
# #   )
#         
# test1 <- sgRNA_fixed_ensembl %>% 
#   filter(
#     str_equal(ensembl_gene_id, new_id) == FALSE
#   ) %>% 
#   dplyr::select(gene, ensembl_gene_id, new_id)

### TEST CODE

# get protospacer 1
# protospacer 1 is followed by constant region 3, based of off pJR89 and pJR101
# construct:
# mU6-protospacer_1-CR3
# hU6-protospacer_2-CR1



feature_ref_df_a <- sgRNA_raw_csv %>% 
  select(unique_sg_rna_pair_id, gene, ensembl_gene_id, sg_id_a, targeting_sequence_a, transcript) %>% 
  mutate(
    # add CR3 as pattern
    pattern = "GTTTCAGAGCTAAGCACAAG"
  )

sg_a_intermediate <- feature_ref_df_a %>% 
  mutate(
    id = paste0(sg_id_a, "-1"), # -1 designates protospacer a
    name = paste0(transcript,"-1"), 
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
    targeting_sequence_b,
    transcript
  ) %>% 
  mutate(
    pattern = "GTTTAAGAGCTAAGCTGGAA"
  )

sg_b_intermediate <- feature_ref_df_b %>% 
  mutate(
    id = paste0(sg_id_b, "-2"), # -1 designates protospacer a
    name = paste0(transcript,"-2"), 
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

feature_ref_csv <- rbind(sg_a_intermediate, sg_b_intermediate) %>% 
  # channge non-targeting values to "Non-Targeting
  mutate(
    across(
      !id,
      \(x) x = case_when(
        str_detect(x, "^non-targeting") ~ "Non-Targeting",
        TRUE ~ x
      )
    )
  ) %>%
  # drop guide RNAs with no ensembl ID
  filter(
    !is.na(target_gene_id)
  ) %>% 
  # replace commas with ensembl ids in feature id with -
  mutate(
    id = str_replace_all(
      id, ",", "-"
    )
  )
  
  
  
  
  
  ### OBSOLETE CODE
  # fetch ensembl id for genes that don't currently have them
  # this code is not good but i cba to figure out something better
  
  # join table with biomaRt ensembl ids with gene name as key
  # left_join(
  #   ensembl_ids,
  #   by = join_by("target_gene_name" == "external_gene_name"),
  #   relationship = "many-to-many"
  # ) %>% 
  # mutate(
  #   target_gene_id = case_when(
  #     # fetch ensembl id for genes that don't currently have them
  #     is.na(target_gene_id) ~ ensembl_gene_id,
  #     TRUE ~ target_gene_id
  #   )
  # ) %>%
  # dplyr::select(-ensembl_gene_id) %>%
  # mutate(
  #   # hack fix to fix 2 last genes that dont get an id
  #   target_gene_id = case_when(
  #     target_gene_name == "AHSA2" ~ "ENSG00000173209",
  #     target_gene_name == "C22orf46" ~ "ENSG00000184208",
  #     TRUE ~ target_gene_id
  #   )
  # ) %>%
  # hack fix to fix left join duplicating rows, probably not the best way
  # but hopefully the ensembl ids point to the same gene symbol even if they change for some reason
  # distinct(
  #   pattern, sequence, .keep_all = TRUE
  # )
  ### END OF OBSOLETE CODE


feature_ref_csv %>% 
  group_by(
    pattern, sequence
  ) %>% 
  mutate(
    n = n()
  ) %>% 
  filter(n > 1)

feature_ref_csv %>% 
  filter(
    if_any(
      everything(),
      \(x) is.na(x)
    )
  )

write_csv(
  feature_ref_csv,
  paste0(
    format(Sys.time(), "%Y%m%d"),
    "_dualguide_sgRNA_feature_ref.csv"
  )
)
