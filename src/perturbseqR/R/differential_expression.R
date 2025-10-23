#' Differential Expression and Feature Selection Utilities
#'
#' Functions for identifying differentially expressed genes and selecting informative
#' features in single-cell datasets.
#'
#' @name differential_expression
NULL

#' Identify overdispersed (noisy) genes
#'
#' Finds genes whose coefficient of variation (CV) exceeds the baseline mean–CV relationship.
#'
#' @param pop A CellPopulation-like object with gene statistics in `pop$genes`.
#' @param noisy_threshold Top quantile of excess CV (default 0.05 for top 5%).
#' @param mean_threshold Minimum mean expression to consider.
#' @param exclude Genes to exclude (vector of names or IDs).
#' @param resolution Number of bins for mean–CV smoothing.
#' @return Character vector of overdispersed gene IDs.
#' @export
find_noisy_genes <- function(pop, noisy_threshold = 0.05, mean_threshold = 0.05,
                             exclude = NULL, resolution = 1000) {
  if (!is.null(exclude)) {
    exclude <- pop$gene_ids(exclude)
  }
  
  thresholded <- subset(pop$genes, mean > mean_threshold)[order(pop$genes$mean, decreasing = TRUE), ]
  gene_means <- thresholded$mean
  gene_cvs <- thresholded$cv
  
  filt_cv <- stats::filter(gene_cvs, rep(1/15, 15))  # smoother than medfilt
  n <- length(gene_means)
  idx <- unique(round(seq(1, n, length.out = resolution)))
  interp_fun <- stats::approxfun(gene_means[idx], filt_cv[idx], rule = 2)
  predicted_cv <- interp_fun(gene_means)
  
  excess_cv <- gene_cvs / predicted_cv
  thresholded$excess_cv <- excess_cv
  
  cv_cutoff <- stats::quantile(excess_cv, probs = 1 - noisy_threshold, na.rm = TRUE)
  noisy <- rownames(subset(thresholded, excess_cv > cv_cutoff & !(rownames(thresholded) %in% exclude)))
  
  message(length(noisy), " variable genes found (", sum(rownames(thresholded) %in% exclude), " excluded)")
  
  plot(1 / sqrt(thresholded$mean), thresholded$cv, pch = 16, cex = 0.4, col = rgb(0,0,0,0.3))
  lines(1 / sqrt(thresholded$mean), predicted_cv, col = "gray")
  points(1 / sqrt(pop$genes[noisy, "mean"]), pop$genes[noisy, "cv"], pch = 16, col = "red")
  
  return(noisy)
}

#' Kolmogorov–Smirnov differential expression
#'
#' Performs two-sample KS tests comparing each subpopulation against a control population.
#'
#' @param pop A CellPopulation-like object.
#' @param key Metadata column defining subpopulations.
#' @param control_cells Query selecting control cells.
#' @param genes Optional subset of genes.
#' @param cells Optional subset of cells.
#' @param normalized Use normalized expression matrix.
#' @param n_jobs Number of cores.
#' @param alpha Significance threshold for multiple testing correction.
#' @param multi_method Method for multiple testing correction (see \code{p.adjust.methods}).
#' @return List of data frames: \code{stat}, \code{p}, and \code{adj_p}.
#' @export
ks_de <- function(pop, key, control_cells, genes = NULL, cells = NULL, normalized = FALSE,
                  n_jobs = 1, alpha = 0.05, multi_method = "BY", ...) {
  ctrl <- pop$where(cells = control_cells, genes = genes, normalized = normalized, ...)
  message(nrow(ctrl), " control cells")
  
  subpops <- pop$groupby(key, cells = cells, genes = genes, normalized = normalized, ...)
  
  worker <- function(name, subpop) {
    ks_vals <- sapply(colnames(subpop), function(g) {
      test <- suppressWarnings(stats::ks.test(subpop[[g]], ctrl[[g]]))
      c(stat = test$statistic, p = test$p.value)
    })
    data.frame(t(ks_vals))
  }
  
  out <- parallel::mclapply(names(subpops), function(n) worker(n, subpops[[n]]), mc.cores = n_jobs)
  stats_df <- do.call(cbind, lapply(out, `[`, "stat"))
  p_df <- do.call(cbind, lapply(out, `[`, "p"))
  
  adj_p <- apply(p_df, 2, p.adjust, method = multi_method)
  list(stat = stats_df, p = p_df, adj_p = adj_p)
}

#' Anderson–Darling differential expression
#'
#' @inherit ks_de
#' @export
ad_de <- function(pop, key, control_cells, genes = NULL, cells = NULL, normalized = FALSE,
                  n_jobs = 1, alpha = 0.05, multi_method = "BY", ...) {
  if (!requireNamespace("goftest", quietly = TRUE))
    stop("Package 'goftest' is required for Anderson–Darling tests.")
  
  ctrl <- pop$where(cells = control_cells, genes = genes, normalized = normalized, ...)
  message(nrow(ctrl), " control cells")
  
  subpops <- pop$groupby(key, cells = cells, genes = genes, normalized = normalized, ...)
  
  worker <- function(name, subpop) {
    ad_vals <- sapply(colnames(subpop), function(g) {
      test <- suppressWarnings(goftest::ad.test(subpop[[g]], ctrl[[g]]))
      c(stat = test$statistic, p = test$p.value)
    })
    data.frame(t(ad_vals))
  }
  
  out <- parallel::mclapply(names(subpops), function(n) worker(n, subpops[[n]]), mc.cores = n_jobs)
  stats_df <- do.call(cbind, lapply(out, `[`, "stat"))
  p_df <- do.call(cbind, lapply(out, `[`, "p"))
  
  adj_p <- apply(p_df, 2, p.adjust, method = multi_method)
  list(stat = stats_df, p = p_df, adj_p = adj_p)
}

#' Random forest feature selector
#'
#' Uses a random forest to rank genes by their ability to discriminate subpopulations.
#'
#' @param pop A CellPopulation-like object.
#' @param key Metadata column defining classes.
#' @param num_features Restrict to top N features.
#' @param genes, cells Gene/cell subsets.
#' @param normalized Whether to use normalized data.
#' @param ignore Genes to exclude.
#' @param n_jobs Parallel threads.
#' @param random_state Optional seed.
#' @return An object of class \code{TreeSelectorResult}.
#' @export
tree_selector <- function(pop, key, num_features = NULL, genes = NULL, cells = NULL,
                          normalized = TRUE, ignore = NULL, n_jobs = 1, random_state = NULL, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' required for tree_selector().")
  
  matrix <- pop$where(cells = cells, genes = genes, normalized = normalized, ...)
  y <- as.factor(pop$cells[rownames(matrix), key])
  
  if (!is.null(ignore)) matrix <- matrix[, !(colnames(matrix) %in% ignore), drop = FALSE]
  
  n_trees <- ceiling(100 * sqrt(ncol(matrix)))
  rf <- randomForest::randomForest(x = matrix, y = y, ntree = n_trees, importance = TRUE, do.trace = TRUE)
  
  importances <- randomForest::importance(rf, type = 1)
  ranked <- sort(importances[, 1], decreasing = TRUE)
  if (!is.null(num_features)) ranked <- head(ranked, num_features)
  
  acc <- mean(rf$predicted == y)
  
  TreeSelectorResult$new(
    classifier = rf,
    selected_genes = names(ranked),
    importances = ranked,
    acc = acc,
    report = capture.output(print(rf)),
    categories = levels(y)
  )
}

#' @title TreeSelectorResult
#' @description R6 class wrapping results from random forest feature selection.
#' @export
TreeSelectorResult <- R6::R6Class("TreeSelectorResult",
                                  public = list(
                                    classifier = NULL,
                                    selected_genes = NULL,
                                    importances = NULL,
                                    acc = NULL,
                                    report = NULL,
                                    categories = NULL,
                                    
                                    initialize = function(classifier, selected_genes, importances, acc, report, categories) {
                                      self$classifier <- classifier
                                      self$selected_genes <- selected_genes
                                      self$importances <- importances
                                      self$acc <- acc
                                      self$report <- paste(report, collapse = "\n")
                                      self$categories <- categories
                                    },
                                    
                                    print = function(...) {
                                      cat(length(self$selected_genes), "differentially expressed features\n")
                                      cat("Feature prediction accuracy:", round(self$acc, 3), "\n\n")
                                      cat(self$report, "\n")
                                    },
                                    
                                    predict = function(matrix) {
                                      x <- matrix[, self$selected_genes, drop = FALSE]
                                      preds <- predict(self$classifier, newdata = x)
                                      factor(preds, levels = self$categories)
                                    },
                                    
                                    score = function(matrix, categories) {
                                      preds <- self$predict(matrix)
                                      mean(preds == categories)
                                    }
                                  )
)