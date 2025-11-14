#' Dimensionality Reduction Transformers
#'
#' PCA, ICA, t-SNE, and UMAP wrappers for single-cell data.
#' @name transformers
NULL

#' PCA Reducer
#' @export
PCAReducer <- R6::R6Class(
  "PCAReducer",
  public = list(
    n_components = NULL,
    whiten = NULL,
    pca_model = NULL,
    loadings_ = NULL,
    reduced_ = NULL,
    explained_variance_ratio_ = NULL,

    initialize = function(n_components = 10, whiten = FALSE) {
      self$n_components <- n_components
      self$whiten <- whiten
    },

    fit = function(X, y = NULL) {
      X <- as.matrix(X)

      pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
      self$pca_model <- pca

      Z <- pca$x[, seq_len(self$n_components), drop = FALSE]
      comps <- pca$rotation[, seq_len(self$n_components), drop = FALSE]

      self$reduced_ <- prepare_dataframe(Z, rownames = rownames(X), prefix = "PCA")
      self$loadings_ <- prepare_dataframe(comps, rownames = colnames(X), prefix = "PCA")

      total_var <- sum(pca$sdev^2)
      selected_var <- pca$sdev[1:self$n_components]^2
      self$explained_variance_ratio_ <- selected_var / total_var
    },

    fit_transform = function(X, y = NULL) {
      self$fit(X)
      self$reduced_
    },

    transform = function(X, y = NULL) {
      X <- as.matrix(X)
      Z <- scale(X, center = self$pca_model$center, scale = FALSE) %*% self$pca_model$rotation
      Z <- Z[, seq_len(self$n_components), drop = FALSE]
      prepare_dataframe(Z, rownames = rownames(X), prefix = "PCA")
    }
  )
)

#' ICA Reducer
#' @export
ICAReducer <- R6::R6Class(
  "ICAReducer",
  public = list(
    n_components = NULL,
    sort_components = NULL,
    reduced_ = NULL,
    mixing_matrix_ = NULL,
    ica_model = NULL,

    initialize = function(n_components = 10, sort_components = TRUE) {
      self$n_components <- n_components
      self$sort_components <- sort_components
    },

    fit = function(X, y = NULL) {
      X <- as.matrix(X)

      ica <- fastICA::fastICA(X, n.comp = self$n_components)
      self$ica_model <- ica

      Z <- ica$S
      mixing <- ica$A

      if (self$sort_components) {
        energy <- apply(mixing, 2, function(x) sqrt(sum(x^2)))
        order <- order(energy, decreasing = TRUE)
        mixing <- mixing[, order, drop = FALSE]
        Z <- Z[, order, drop = FALSE]
      }

      self$reduced_ <- prepare_dataframe(Z, rownames = rownames(X), prefix = "ICA")
      self$mixing_matrix_ <- prepare_dataframe(mixing, rownames = colnames(X), prefix = "ICA")
    },

    fit_transform = function(X, y = NULL) {
      self$fit(X)
      self$reduced_
    },

    transform = function(X, y = NULL) {
      stop("ICA transform not implemented to match original Python design.")
    }
  )
)

#' PCA + TSNE reducer
#' @export
PCATSNEReducer <- R6::R6Class(
  "PCATSNEReducer",
  public = list(
    n_components = NULL,
    use_pca = NULL,
    pca_model = NULL,
    pca_matrix_ = NULL,
    reduced_ = NULL,

    #' @description
    #' PCA reducer mimicking PCATSNEReducer but *without* any t-SNE.
    initialize = function(n_components = 10, use_pca = TRUE) {
      self$n_components <- n_components
      self$use_pca <- use_pca
    },

    #' @description
    #' Fit PCA model
    fit = function(X, y = NULL) {
      X <- as.matrix(X)

      if (self$use_pca) {
        message("Performing PCA...")
        pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
        self$pca_model <- pca

        Y <- pca$x[, seq_len(self$n_components), drop = FALSE]
        self$pca_matrix_ <- Y
      } else {
        # passthrough
        Y <- X[, seq_len(self$n_components), drop = FALSE]
        self$pca_matrix_ <- Y
        self$pca_model <- NULL
      }

      # reduced_ matches the structure of TSNE reducer's output
      self$reduced_ <- prepare_dataframe(
        self$pca_matrix_,
        rownames = rownames(X),
        prefix = "PCA"
      )
    },

    #' @description
    #' Fit and return reduced coordinates
    fit_transform = function(X, y = NULL) {
      self$fit(X)
      self$reduced_
    },

    #' @description
    #' Apply previously computed PCA rotation
    transform = function(X, y = NULL) {
      if (is.null(self$pca_model)) {
        stop("transform() requires PCA to have been fit with use_pca=TRUE.")
      }
      X <- as.matrix(X)

      Z <- scale(X,
                 center = self$pca_model$center,
                 scale = FALSE) %*%
           self$pca_model$rotation[, seq_len(self$n_components), drop = FALSE]

      prepare_dataframe(Z, rownames = rownames(X), prefix = "PCA")
    }
  )
)

#' PCA + UMAP reducer
#' @export
PCAUMAPReducer <- R6::R6Class(
  "PCAUMAPReducer",
  public = list(
    n_components = NULL,
    n_neighbors = NULL,
    metric = NULL,
    use_pca = NULL,
    pca_model = NULL,
    reduced_ = NULL,
    pca_matrix_ = NULL,
    graph_ = NULL,

    initialize = function(n_components = 10,
                          n_neighbors = 10,
                          metric = "euclidean",
                          use_pca = TRUE) {
      self$n_components <- n_components
      self$n_neighbors <- n_neighbors
      self$metric <- metric
      self$use_pca <- use_pca
    },

    fit = function(X, y = NULL) {
      X <- as.matrix(X)

      if (self$use_pca) {
        message("Performing PCA...")
        pca <- stats::prcomp(X)
        Y <- pca$x[, seq_len(self$n_components), drop = FALSE]
        self$pca_model <- pca
        self$pca_matrix_ <- Y
      } else {
        Y <- X
      }

      message("Performing UMAP...")
      Z <- uwot::umap(Y,
                      n_neighbors = self$n_neighbors,
                      metric = self$metric,
                      verbose = FALSE)

      self$reduced_ <- prepare_dataframe(Z, rownames = rownames(X), prefix = "UMAP")
    },

    fit_transform = function(X, y = NULL) {
      self$fit(X)
      self$reduced_
    },

    transform = function(X, y = NULL) {
      stop("UMAP transform not implemented to match Python behavior.")
    }
  )
)