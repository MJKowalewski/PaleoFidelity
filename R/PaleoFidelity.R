#' PaleoFidelity: Measuring paleontological fidelity (live-dead community agreement)
#'
#' PaleoFidelity package provides common measures of paleontological fidelity
#' used to assess live-dead community agreement. These measures encompass
#' multiple aspects of fidelity, including faunal composition and alpha diversity.
#'
#' Plot functions to visualize fidelity patterns are also provided.
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
#' FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, iter=49)
#'
#'
#' @importFrom vegan vegdist rrarefy
NULL
