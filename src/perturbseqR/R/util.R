#' Utility Functions
#'
#' General mathematical and matrix helper functions.
#' @name util
NULL

#' Non-zero upper triangle
#'
#' Return the non-zero part of the upper triangle of a matrix as a vector.
#'
#' @param M A numeric matrix.
#' @param k Diagonal offset (default = 1). 1 excludes the main diagonal.
#' @return A numeric vector of the nonzero upper-triangular elements.
#' @export
nzflat <- function(M, k = 1) {
  stopifnot(is.matrix(M))
  keep <- upper.tri(M, diag = FALSE)
  if (k > 1) keep <- upper.tri(M, diag = FALSE, k = k)
  vals <- M[keep]
  vals[vals != 0]
}

#' Upper triangle as vector
#'
#' Return the upper triangular part of a matrix in stacked (vector) format.
#'
#' @param M A numeric matrix.
#' @param k Diagonal offset (default = 1).
#' @return A numeric vector corresponding to the upper-triangular portion of M.
#' @export
upper_triangle <- function(M, k = 1) {
  stopifnot(is.matrix(M))
  keep <- upper.tri(M, diag = FALSE)
  if (k > 1) keep <- upper.tri(M, diag = FALSE, k = k)
  as.vector(M[keep])
}

#' Convert categorical columns to strings
#'
#' Replace factor (categorical) columns in a data frame with their character equivalents.
#' This is useful before writing to formats like HDF5 that do not support R factors.
#'
#' @param df A data frame.
#' @return A data frame with all factor columns converted to character type.
#' @export
.strip_cat_cols <- function(df) {
  stopifnot(is.data.frame(df))
  cat_cols <- names(df)[sapply(df, is.factor)]
  if (length(cat_cols) > 0) {
    message("! Converting categorical columns to string...")
    for (col in cat_cols) {
      df[[col]] <- as.character(df[[col]])
    }
  }
  df
}

#' Gini coefficient
#'
#' Compute the Gini coefficient of inequality for a numeric vector.
#'
#' @param x A numeric vector (non-negative values only).
#' @return A single numeric value representing the Gini coefficient.
#' @examples
#' gini(c(1, 2, 3, 4, 5))
#' @export
gini <- function(x) {
  stopifnot(is.numeric(x))
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  x <- x - min(x)
  x <- x + 1e-7
  x <- sort(x)
  n <- length(x)
  index <- seq_len(n)
  (sum((2 * index - n - 1) * x)) / (n * sum(x))
}