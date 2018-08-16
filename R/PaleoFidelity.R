#' PaleoFidelity: Statistical assessment of paleontological fidelity (live-dead community agreement)
#'
#' PaleoFidelity package provides common measures of paleontological fidelity (live-dead community
#' agreement), including compositional agreement, alpha diversity, and beta diversity.
#'
#' PaleoFidelity offers permutation tests for one-sample and multi-sample tests, and
#' bootstrapped confidence intervals for fidelity statistics.
#'
#' Several plot functions are included to visualize outcomes of fidelity analyses.
#'
#' @author Michal Kowalewski \email{kowalewski@@ufl.edu}
#' @name PaleoFidelity-package
#' @aliases PaleoFidelity-package PaleoFidelity
#' @docType package
#'
#' @examples
#' # Version and citation
#' packageVersion('PaleoFidelity')
#' citation('PaleoFidelity')
#'
#' # Check if data are compliant and generate basic data summary
#' \dontrun{require(vegan)
#' data(dune)
#' FidelitySummary(as.matrix(dune), as.matrix(dune[sample(1:nrow(dune)),]), report=T)
#'
#' # Compute measures of compositional fidelity
#' out1 <- FidelityEst(as.matrix(dune), as.matrix(dune[sample(1:nrow(dune)),]))
#' plot(out1)}
NULL
