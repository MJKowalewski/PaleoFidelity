#' PaleoFidelity: Measuring paleontological fidelity (live-dead community agreement)
#'
#' PaleoFidelity package provides common measures of paleontological fidelity
#' used to assess live-dead community agreement. These measures encompass
#' multiple aspects of fidelity, including faunal composition and alpha diversity.
#' Plot functions to visualize fidelity patterns are also provided.
#'
#' The PaleoFidelity functions are designed to evaluate and
#' visualize congruence (fidelity) between living communities
#' and sympatric death assemblages (the incipient fossil
#' record). The intellectual motivation for live-dead fidelity analyses
#' is the assessments of the informative value of the fossil record
#' ('taphonomy') and the evaluation of recent human-induced ecosystem
#' changes reflected in live-dead disagreements ('conservation
#' paleobiology').
#'
#' The current version of the package provides functions for
#' assessing 'compositional fidelity' (the agreement in community
#' composition in terms of correlation and similarity measures) and
#' 'diversity fidelity' (offsets in species richness and evenness).
#' The functions include sample standardizations, resampling models,
#' and by-group analyses when grouping factors are provided. Plotting
#' functions designed to visualize compositional and diversity fidelity
#' are also included.
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
#' FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat,
#' n.filters=30, report=TRUE)
#'
#' # Compositional fidelity analysis with sites grouped by habitat
#' data(FidData)
#' out1 <- FidelityEst(live = FidData$live[6:9,], dead = FidData$dead[6:9,],
#'                     gp = FidData$habitat[6:9], cor.measure='spearman',
#'                     sim.measure='bray', n.filters=30, iter=99, rm.zero=FALSE, tfsd='none')
#' SJPlot(out1, gpcol=c('forestgreen', 'coral3'))
#'
#' # Diversity fidelity analysis with sites grouped by habitat
#' my.fid <- FidelityDiv(FidData$live[6:9,], FidData$dead[6:9,],
#'                       FidData$habitat[6:9], n.filters=20,
#'                        iter=100, CI=0.95)
#' # site-level estimates of Delta S with 95% CIs and p values
#' my.fid$x
#' # p values for means of groups
#' my.fid$p.gps
#' AlphaPlot(my.fid, col.gp=c('forestgreen', 'coral1'), bgpt='beige')
#'
#' @importFrom vegan vegdist rrarefy
NULL
