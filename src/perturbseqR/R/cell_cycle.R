#' Cell cycle analysis utilities
#'
#' Functions to compute and visualize cell cycle phase and progression from single-cell
#' expression data, adapted from the original Perturb-seq Python implementation.
#'
#' @name cell_cycle
NULL

#' Group correlation within a gene set
#'
#' Measures the correlation of all genes within a list to the mean expression
#' of all genes in that list.
#'
#' @param population A CellPopulation-like object (must support `where()` returning a matrix or data frame of expression values).
#' @param gene_list A character vector of gene names.
#' @return A numeric named vector of correlation coefficients for each gene.
#' @export
group_corr <- function(population, gene_list) {
  expression_matrix <- population$where(genes = gene_list)
  total <- rowMeans(expression_matrix, na.rm = TRUE)
  corrs <- apply(expression_matrix, 2, function(col) cor(col, total, use = "pairwise.complete.obs"))
  corrs
}

#' Refine a gene list by correlation threshold
#'
#' Removes genes that do not correlate strongly with the mean expression of the group.
#'
#' @param population A CellPopulation-like object.
#' @param gene_list Character vector of gene names.
#' @param threshold Minimum correlation coefficient to retain a gene.
#' @param return_corrs Logical; whether to return correlations as well as gene names.
#' @return A character vector of refined genes, or a data frame if `return_corrs = TRUE`.
#' @export
refine_gene_list <- function(population, gene_list, threshold = 0.2, return_corrs = FALSE) {
  corrs <- group_corr(population, gene_list)
  keep <- corrs >= threshold
  if (return_corrs) {
    data.frame(gene = names(corrs), correlation = corrs)[keep, ]
  } else {
    names(corrs)[keep]
  }
}

#' Compute a cell-cycle group score
#'
#' Z-normalized log-expression score of a given gene set across all cells.
#'
#' @param population A CellPopulation-like object.
#' @param gene_list Character vector of gene names.
#' @return A numeric vector of z-scored gene set expression per cell.
#' @export
group_score <- function(population, gene_list) {
  expr <- population$where(genes = gene_list)
  log_expr <- log2(expr + 1)
  scores <- rowSums(log_expr, na.rm = TRUE)
  scale(scores)
}

#' Compute scores for multiple gene sets
#'
#' @param population A CellPopulation-like object.
#' @param gene_lists Named list of character vectors of gene names.
#' @return A data frame of z-scored expression for each gene set.
#' @export
batch_group_score <- function(population, gene_lists) {
  out <- lapply(gene_lists, function(g) group_score(population, g))
  as.data.frame(out)
}

#' Get canonical cell-cycle gene sets
#'
#' Returns a named list of known marker genes for cell-cycle phases.
#'
#' @param pop Optional population for refinement.
#' @param refine Logical; refine based on correlation.
#' @param threshold Minimum correlation to keep gene if refining.
#' @return A named list of character vectors.
#' @export
get_cell_phase_genes <- function(pop = NULL, refine = FALSE, threshold = 0.3) {
  phase_genes <- list(
    "G1-S" = c("ARGLU1", "BRD7", "CDC6", "CLSPN", "ESD", "GINS2",
               "GMNN", "LUC7L3", "MCM5", "MCM6", "NASP", "PCNA",
               "PNN", "SLBP", "SRSF7", "SSR3", "ZRANB2"),
    "S" = c("ASF1B", "CALM2", "CDC45", "CDCA5", "CENPM", "DHFR",
            "EZH2", "FEN1", "HIST1H2AC", "HIST1H4C", "NEAT1", "PKMYT1",
            "PRIM1", "RFC2", "RPA2", "RRM2", "RSRC2", "SRSF5", "SVIP",
            "TOP2A", "TYMS", "UBE2T", "ZWINT"),
    "G2-M" = c("AURKB", "BUB3", "CCNA2", "CCNF", "CDCA2", "CDCA3",
               "CDCA8", "CDK1", "CKAP2", "DCAF7", "HMGB2", "HN1",
               "KIF5B", "KIF20B", "KIF22", "KIF23", "KIFC1", "KPNA2",
               "LBR", "MAD2L1", "MALAT1", "MND1", "NDC80", "NUCKS1",
               "NUSAP1", "PIF1", "PSMD11", "PSRC1", "SMC4", "TIMP1",
               "TMEM99", "TOP2A", "TUBB", "TUBB4B", "VPS25"),
    "M" = c("ANP32B", "ANP32E", "ARL6IP1", "AURKA", "BIRC5", "BUB1",
            "CCNA2", "CCNB2", "CDC20", "CDC27", "CDC42EP1", "CDCA3",
            "CENPA", "CENPE", "CENPF", "CKAP2", "CKAP5", "CKS1B",
            "CKS2", "DEPDC1", "DLGAP5", "DNAJA1", "DNAJB1", "GRK6",
            "GTSE1", "HMG20B", "HMGB3", "HMMR", "HN1", "HSPA8",
            "KIF2C", "KIF5B", "KIF20B", "LBR", "MKI67", "MZT1",
            "NUF2", "NUSAP1", "PBK", "PLK1", "PRR11", "PSMG3", "PWP1",
            "RAD51C", "RBM8A", "RNF126", "RNPS1", "RRP1", "SFPQ",
            "SGOL2", "SMARCB1", "SRSF3", "TACC3", "THRAP3", "TPX2",
            "TUBB4B", "UBE2D3", "USP16", "WIBG", "YWHAH", "ZNF207"),
    "M-G1" = c("AMD1", "ANP32E", "CBX3", "CDC42", "CNIH4", "CWC15",
               "DKC1", "DNAJB6", "DYNLL1", "EIF4E", "FXR1", "GRPEL1",
               "GSPT1", "HMG20B", "HSPA8", "ILF2", "KIF5B", "KPNB1",
               "LARP1", "LYAR", "MORF4L2", "MRPL19", "MRPS2", "MRPS18B",
               "NUCKS1", "PRC1", "PTMS", "PTTG1", "RAN", "RHEB", "RPL13A",
               "SRSF3", "SYNCRIP", "TAF9", "TMEM138", "TOP1", "TROAP",
               "UBE2D3", "ZNF593")
  )
  if (refine && !is.null(pop)) {
    for (phase in names(phase_genes)) {
      phase_genes[[phase]] <- refine_gene_list(pop, phase_genes[[phase]], threshold)
    }
  }
  phase_genes
}

#' Compute cell-cycle phase scores
#'
#' @param pop A CellPopulation-like object.
#' @param gene_list Optional named list of marker genes.
#' @param refine Logical; whether to refine gene sets.
#' @param threshold Minimum correlation coefficient for refinement.
#' @return A data frame of phase scores and assignments.
#' @export
get_cell_phase <- function(pop, gene_list = NULL, refine = TRUE, threshold = 0.3) {
  if (is.null(gene_list))
    gene_list <- get_cell_phase_genes(pop, refine = refine, threshold = threshold)
  
  phase_scores <- batch_group_score(pop, gene_list)
  norm_scores <- scale(phase_scores)
  
  cor_scores <- cor(t(norm_scores))
  phase_list <- names(gene_list)
  
  # Assign phase with maximal correlation
  max_phase <- apply(cor_scores, 1, function(row) phase_list[which.max(row)])
  cell_cycle_scores <- as.data.frame(norm_scores)
  cell_cycle_scores$cell_cycle_phase <- factor(max_phase, levels = phase_list)
  cell_cycle_scores
}

#' Add cell-cycle scores to a population
#'
#' @param pop A CellPopulation-like object.
#' @param gene_list Optional named list of marker genes.
#' @return The population with appended cell cycle annotations.
#' @export
add_cell_cycle_scores <- function(pop, gene_list = NULL) {
  scores <- get_cell_phase(pop, gene_list = gene_list)
  pop$add_property(cells = scores)
  pop$cells$cell_cycle_phase <- factor(pop$cells$cell_cycle_phase,
                                       levels = names(get_cell_phase_genes()))
  pop
}

#' Plot a cell-cycle position heatmap
#'
#' Visualize ordered cells by their phase and progress through the cell cycle.
#'
#' @param pop A CellPopulation-like object.
#' @param cells Optional subset of cells.
#' @param ... Additional arguments passed to internal query.
#' @export
cell_cycle_position_heatmap <- function(pop, cells = NULL, ...) {
  if (is.null(cells)) {
    df <- pop$cells[, c("G1-S", "S", "G2-M", "M", "M-G1", "cell_cycle_phase", "cell_cycle_progress")]
  } else {
    celllist <- pop$where(cells = cells, ...)$index
    df <- pop$cells[celllist, c("G1-S", "S", "G2-M", "M", "M-G1", "cell_cycle_phase", "cell_cycle_progress")]
  }
  df <- df[order(df$cell_cycle_phase, -df$cell_cycle_progress), ]
  mat <- t(as.matrix(df[, 1:5]))
  heatmap(mat, Rowv = NA, Colv = NA, scale = "none", col = hcl.colors(50, "Blu-Yl"))
}
