#' PerturbseqR initialization
#'
#' Startup configuration and metadata for the PerturbseqR package.
#' @noRd

.PERTURBSEQ_VERSION <- "0.1"

.onLoad <- function(libname, pkgname) {
  packageStartupMessage(paste("PerturbseqR version", .PERTURBSEQ_VERSION, "loaded."))
}