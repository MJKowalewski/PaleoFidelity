#' FidelityEst (computes live-dead fidelity measures)
#'
#' FidelityEst computes common live-dead fidelity measures for two datasets (live and dead)
#'
#' @param live A matrix with counts of live-collected specimens (rows=samples, columns=species/taxa)
#'
#' @param dead A matrix with counts of dead specimens (rows=samples, columns=species/taxa)
#'    (dimensions must match)
#'
#' @param gp An optional factor, with two or more levels, defining sample groups
#'    (the length of gp must equal number of rows in live and dead)
#'
#' @param cor.measure Correlation coefficient to measure pairwise live-dead correlations (default='Spearman')
#'
#' @param dist.measure Measure of similarity to measure live-dead agreement (default='Chao')
#'
#' @return A list containing the following components:
#'   \item{x}{The values of the Spearman coefficients \emph{rho} for each live-dead comparison}
#'
#' @examples
#'
#' library(vegan)
#' data(dune)
#' \dontrun{
#' out1 <- FidelityEst(as.matrix(dune), as.matrix(dune[sample(1:nrow(dune)),]))
#' plot(out1, xlim=c(-1,1), ylim=c(0,1), las=1, xlab=out1$xlab, ylab=out1$ylab)
#' abline(h=0.5, v=0)}


FidelityEst <- function(live, dead, gp=NULL, cor.measure='spearman', dist.measure='chao') {
  FidelitySummary(live, dead, gp)
  x1 <- as.data.frame(t(live))
  x2 <- as.data.frame(t(dead))
  cor.e <- mapply(function(x,y) stats::cor(x,y, method=cor.measure), x1, x2)
  dist.e <- mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), 'chao'), x1, x2)
  out1 <- list(x=cor.e, y=dist.e, xlab=cor.measure, ylab=dist.measure)
  return(out1)
}
