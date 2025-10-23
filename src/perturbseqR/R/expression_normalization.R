#' Expression normalization utilities for Perturb-seq data
#' 
#' Functions for normalizing and transforming single-cell expression matrices
#' within a CellPopulation-like object.
#' 
#' @name expression_normalization
NULL

#' @import dplyr Matrix
#' @importFrom stats median sd
#' @importFrom utils head tail

#' Remove low-expressing genes
#' 
#' Removes genes with mean expression below a threshold to reduce memory usage.
#' @param pop A CellPopulation-like object (with `matrix` and `genes` slots)
#' @param threshold Minimum mean expression to retain
#' @export
strip_low_expression <- function(pop, threshold = 0) {
  retain <- rownames(pop$genes[pop$genes$mean > threshold, ])
  remove <- rownames(pop$genes[pop$genes$mean <= threshold, ])
  
  if (length(remove) == 0) {
    message("No genes have expression below threshold.")
    return(pop)
  }
  
  pop$matrix <- pop$matrix[retain, , drop = FALSE]
  pop$genes[remove, "in_matrix"] <- FALSE
  
  numeric_cols <- names(pop$genes)[sapply(pop$genes, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("gene_name", "in_matrix"))
  
  for (col in numeric_cols) {
    pop$genes[remove, col] <- NA_real_
  }
  
  gc()
  pop
}

#' Normalize all cells to equal UMI counts
#'
#' Scales each cell’s expression vector to the same median UMI count.
#' @param matrix Expression matrix (cells × genes)
#' @param median_umi_count Optional fixed UMI normalization target
#' @export
equalize_UMI_counts <- function(matrix, median_umi_count = NULL) {
  reads_per_bc <- rowSums(matrix)
  if (is.null(median_umi_count)) {
    median_reads_per_bc <- median(reads_per_bc)
  } else {
    median_reads_per_bc <- median_umi_count
  }
  scaling_factors <- median_reads_per_bc / reads_per_bc
  m <- sweep(matrix, 1, scaling_factors, "*")
  
  if (mean(median_reads_per_bc) < 5000) {
    message("Scaling with a small number of reads. Are you sure this is what you want?")
  }
  m
}

#' Log-normalize expression
#' 
#' Applies log2(X + pseudocount) normalization, optionally after UMI equalization.
#' @param pop CellPopulation-like object
#' @param scale_by_total Whether to equalize UMI counts
#' @param pseudocount Offset for zero values
#' @return Data frame of log-normalized expression
#' @export
log_normalize_expression <- function(pop, scale_by_total = TRUE, pseudocount = 1) {
  matrix <- pop$matrix
  if (scale_by_total) {
    m <- equalize_UMI_counts(matrix)
  } else {
    m <- as.matrix(matrix)
  }
  log2(m + pseudocount)
}

#' Z-normalize expression
#'
#' Normalizes each gene across cells by Z-score.
#' @param pop CellPopulation-like object
#' @param scale_by_total Whether to equalize UMI counts first
#' @return Data frame of Z-normalized expression
#' @export
z_normalize_expression <- function(pop, scale_by_total = TRUE) {
  matrix <- pop$matrix
  if (scale_by_total) {
    m <- equalize_UMI_counts(matrix)
  } else {
    m <- as.matrix(matrix)
  }
  m_out <- scale(m, center = TRUE, scale = TRUE)
  as.data.frame(m_out)
}

#' Normalize expression to control population
#'
#' Scales a test matrix relative to a control matrix (e.g., DMSO cells).
#' @param matrix Expression matrix
#' @param control_matrix Control expression matrix
#' @param scale_by_total Whether to equalize UMIs first
#' @param median_umi_count Optional fixed UMI normalization target
#' @export
normalize_matrix_to_control <- function(matrix, control_matrix, 
                                        scale_by_total = TRUE, median_umi_count = NULL) {
  if (inherits(matrix, "dgCMatrix")) matrix <- as.matrix(matrix)
  if (inherits(control_matrix, "dgCMatrix")) control_matrix <- as.matrix(control_matrix)
  
  if (scale_by_total) {
    message("     Determining scale factors...")
    reads_per_bc <- rowSums(matrix)
    median_reads_per_bc <- if (is.null(median_umi_count)) median(reads_per_bc) else median_umi_count
    scaling_factors <- median_reads_per_bc / reads_per_bc
    
    message("     Normalizing matrix to median")
    m <- sweep(matrix, 1, scaling_factors, "*")
    
    if (mean(median_reads_per_bc) < 5000) {
      message("Scaling with a small number of reads. Are you sure this is what you want?")
    }
    
    control_reads_per_bc <- rowSums(control_matrix)
    control_scaling_factors <- median_reads_per_bc / control_reads_per_bc
    message("     Normalizing control matrix to median")
    c_m <- sweep(control_matrix, 1, control_scaling_factors, "*")
  } else {
    m <- as.matrix(matrix)
    c_m <- as.matrix(control_matrix)
  }
  
  control_mean <- colMeans(c_m)
  control_sd <- apply(c_m, 2, sd)
  
  message("     Scaling matrix to control")
  m_out <- sweep(m, 2, control_mean, "-")
  m_out <- sweep(m_out, 2, control_sd, "/")
  
  message("     Done.")
  as.data.frame(m_out)
}

#' Inherit normalized matrix from parent population
#' 
#' Copies normalized matrix values from a parent CellPopulation subset.
#' @param pop Target subpopulation
#' @param parent_pop Parent population with normalized_matrix slot
#' @export
inherit_normalized_matrix <- function(pop, parent_pop) {
  pop$normalized_matrix <- parent_pop$normalized_matrix[rownames(pop$matrix), colnames(pop$matrix)]
  pop
}