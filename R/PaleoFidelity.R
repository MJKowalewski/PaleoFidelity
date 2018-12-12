#' PaleoFidelity: Measuring paleontological fidelity (live-dead community agreement)
#'
#' PaleoFidelity package provides common measures of paleontological fidelity
#' used to assess live-dead community agreement. These measures encompass
#' multiple aspects of fidelity, including faunal composition, alpha diversity,
#' and beta diversity.
#'
#' PaleoFidelity package also offers permutation tests for one-sample, two-sample, and
#' multi-sample hypotheses as well as bootstrapped confidence intervals for fidelity statistics.
#'
#' Specialized plot functions are included to visualize fidelity patterns.
#'
#' @author Michal Kowalewski \email{kowalewski@@ufl.edu}
#' @name PaleoFidelity
#' @aliases PaleoFidelity
#' @docType package
#'
#' @examples
#' # Version and citation
#' packageVersion('PaleoFidelity')
#' citation('PaleoFidelity')
#'
#' # Check if data are compliant and generate basic data summary
#' data(FidData)
#' FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, report=TRUE)
#'
#' # Returns estimates of compositional fidelity
#' FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat)
#'
#'
#' @importFrom vegan vegdist rrarefy
NULL
