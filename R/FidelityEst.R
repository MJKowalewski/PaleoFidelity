#' FidelityEst computes compositional measures of live-dead fidelity
#'
#' FidelityEst computes common live-dead fidelity measures by comparing two matching
#' matrices (live and dead) with species/taxon abundance data.
#' In the case of datastes representing more than one sample,
#' the function returns also means. If 'gp' factor is provided to define sample groups
#' means for groups are returned as well.
#'
#'
#' @param live A species/taxon abundance matrix with counts of live-collected specimens
#'    (rows=samples, columns=species/taxa)
#'
#' @param dead A matrix with counts of dead specimens (rows=samples, columns=species/taxa)
#'    (dimensions of 'live' and 'dead' must match)
#'
#' @param gp An optional factor, with two or more levels, defining sample groups
#'    (the length of gp must equal number of rows in live and dead)
#'
#' @param cor.measure defines correlation coefficient (stats function 'vignettecor') used to measure pairwise
#'    live-dead correlations (default='spearman').
#'
#' @param sim.measure A measure of similarity to measure live-dead agreement (default='Chao')
#'
#' @param report Print compliance report from function FidelitySummary
#'
#' @return A list containing the following components:
#'   \item{x}{The values of correlation coefficient for each live-dead comparison}
#'   \item{y}{The values of similarity coefficients for each live-dead comparison}
#'   \item{x.lab}{The name of the correlation coefficient used}
#'   \item{y.lab}{The name of the similarity coefficient used}
#'   \item{mean.x}{The mean correlation coefficient (same as 'x' if live-dead data include only one sample )}
#'   \item{mean.y}{The mean similarity coefficient (same as 'x' if live-dead data include only one sample )}
#'   \item{group.means.x}{Mean correlation coefficients by group (returned if 'gp' factor provided)}
#'   \item{group.means.y}{Mean similarity coefficients by group (returned if 'gp' factor provided)}
#'
#' @examples
#'
#' data(FidData)
#' FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat)
#'
#' @export
#' @importFrom stats sd
#' @importFrom vegan vegdist

FidelityEst <- function(live, dead, gp=NULL, cor.measure='spearman', sim.measure='chao', report=FALSE) {
  FidelitySummary(live, dead, gp, report=report)
  x1 <- as.data.frame(t(live))
  x2 <- as.data.frame(t(dead))
  cor.e <- mapply(function(x,y) stats::cor(x,y, method=cor.measure), x1, x2)
  sim.e <- mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), dist=sim.measure), x1, x2)
  mean.cor <- mean(cor.e)
  mean.sim <- mean(sim.e)
  if(length(gp) == 0)
   out1 <- list(x=cor.e, y=sim.e,
                x.lab=cor.measure, y.lab=sim.measure,
                mean.x=mean.cor, mean.y=mean.sim)
  if(length(gp) > 0)
   mean.cor.gp <- tapply(cor.e, gp, mean)
   mean.sim.gp <- tapply(sim.e, gp, mean)
   out1 <- list(x=cor.e, y=sim.e,
                x.lab=cor.measure, y.lab=sim.measure,
                mean.x=mean.cor, mean.y=mean.sim,
                group.means.x=mean.cor.gp, group.means.y=mean.sim.gp)
  return(out1)
}
