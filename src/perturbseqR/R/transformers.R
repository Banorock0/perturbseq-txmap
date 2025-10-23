#' Dimensionality Reduction Transformers
#'
#' PCA, ICA, t-SNE, and UMAP wrappers for single-cell data.
#' @name transformers
NULL

#' PCA Reducer
#' @export
PCAReducer <- function(expr, n_components = 10) {
  stop("PCAReducer() not yet implemented.")
}

#' ICA Reducer
#' @export
ICAReducer <- function(expr, n_components = 10) {
  stop("ICAReducer() not yet implemented.")
}

#' PCA + t-SNE Reducer
#' @export
PCATSNEReducer <- function(expr, n_components = 2) {
  stop("PCATSNEReducer() not yet implemented.")
}

#' PCA + UMAP Reducer
#' @export
PCAUMAPReducer <- function(expr, n_components = 2) {
  stop("PCAUMAPReducer() not yet implemented.")
}
