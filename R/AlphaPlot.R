#' AlphaPlot (outputs a cross plot of live-dead offsets in alpha diversity and evenness)
#'
#' AlphaPlot function generates a bivariate plot of an offset in alpha diversity (DELTA S)
#' on x axis versus an offset in evenness (DELTA PIE) for all live-dead pairwise comparisons.
#'
#' @details AlphaPlot function uses the output of FidelityDiv function to produce
#' a cross plot visualizing differences in estimates of alpha diversty and evenness between
#' pairs of sympatric live and dead samples. In its default form, an approch proposed by
#' Olszewski and Kidwell (200X) is used. The offset in alpha diversity is measured as the difference
#' between natural logarithms of sample standardized richness log(R-dead) - log(R-live). The offset
#' in evenness is measured as the difference between Hurlbert's PIE (PIE-dead - PIE-live).
#' Confidence bars depitcting user-specified confidence intervals are also plotted.
#' If a grouping factor is provided, symbols and bars are color-coded by levels and
#' group means are plotted in addition to grand mean for all sites.
#'
#'
#' @param x An object (a list) returned by FidelityDiv function
#'
#' @param CI A numerical value (default = 0.5) defining confidence bars for individual sites.
#'  Note: 0.5 - plots bars represnting interuqartile range, 0.95 - plots 95% confidence bars, etc.
#'  Confidence bars are estimated as percentiles of subsampled estimates of Delta S and Delta PIE.
#'
#' @param CImean A numerical value (default = 0.99) defining confidence bars for means of all sites
#'  or groups of sites (if 'gp' factor was provided). Note: 0.5 - plots bars represnting interuqartile
#'  range, 0.95 - plots 95% confidence bars, etc. Confidence bars are estimated as percentiles of
#'  subsampled estimates of Delta S and Delta PIE.
#'
#' @param colpt A character string or number defining color of symbols for individual sites
#' (default = 'gray')
#'
#' @param colci A character string or number defining color of confidence bars for individual
#'  sites (default = 'gray')
#'
#' @param bgpt A character string or number defining color of symbols for individual
#' sites (default = 'white')
#'
#' @param pchpt A number defining symbol type (default = 21)
#' Note: filled symbol numbers (21-25) should be used
#'
#' @param colmean A character string or number defining color of bars representing
#' mean estimates for all sites (default = 'white')
#'
#'
#' @param colpt.gp A character string or number defining colors by group (default = 1:#levels)
#' for symbols for sites (applicable when 'gp' factor is provided). Custom colors
#' can be provided, but must match number of levels in a grouping factor
#'
#' @param colci.gp A character string or number defining colors by group (default = 1:#levels)
#' for confidence bars for sites (applicable when 'gp' factor is provided). Custom colors
#' can be provided, but must match number of levels in a grouping factor
#'
#'
#' @param addlegend Logical (default=TRUE): add legend (activated only when 'gp' factor provided)
#'
#'
#' @return A single bivariate plot produce by plot function
#'
#'
#' @examples
#'
#' AlphaPlot(FidelityDiv(FidData$live, FidData$dead, FidData$habitat, n.filters=30, t.filters=1))
#'
#' @export

AlphaPlot <- function(x, CI=0.5, CImean=0.99, colpt='gray', colci='gray', bgpt='white',
                      pchpt=21, colmean='black', colpt.gp=NULL, colci.gp=NULL,
                      addlegend=TRUE) {
  z <- x$out
  samDS <- rowMeans(z[,1,])
  samDScf <- t(apply(z[,1,], 1, stats::quantile, prob=c((1-CI)/2, 1 - (1-CI)/2)))
  samDP <- rowMeans(z[,2,])
  samDPcf <- t(apply(z[,2,], 1, stats::quantile, prob=c((1-CI)/2, 1 - (1-CI)/2)))
  meanDS <- mean(colMeans(z[,1,]))
  meanDP <- mean(colMeans(z[,2,]))
  DScf <- stats::quantile(colMeans(z[,1,]), prob=c((1-CImean)/2, 1 - (1-CImean)/2))
  DPcf <- stats::quantile(colMeans(z[,2,]), prob=c((1-CImean)/2, 1 - (1-CImean)/2))
  xmax <- max(abs(range(samDScf)))
  ymax <- max(abs(range(samDPcf)))

  graphics::plot(samDS, samDP, xlim=c(-xmax, xmax), ylim=c(-ymax, ymax), type='n', las=1,
       xlab=bquote(Delta['S']), ylab=bquote(Delta['PIE']), cex.lab=1.5)
  graphics::abline(h=0, v=0, lwd=1, col='gray')
  if (length(x$gp)==0) {
   for(i in 1:length(samDS)) {
     graphics::points(cbind(cbind(rep(samDS[i], 2)), cbind(samDPcf[i,])), type='l', col=colci, lwd=0.5, lend=3)
     graphics::points(cbind(cbind(samDScf[i,]), cbind(rep(samDP[i],2))), type='l', col=colci, lwd=0.5, lend=3)
     }
    graphics::points(samDS, samDP, pch=pchpt, col=colpt, bg=bgpt, cex=0.8)
  }
  if (length(x$gp)>0) {
    ifelse(length(colci.gp) > 0, colci.gp <- colci.gp, colci.gp <- 1:(length(levels(x$gp))))
    barcols <- colci.gp[x$gp]
    for(i in 1:length(samDS)) {
      graphics::points(cbind(cbind(rep(samDS[i], 2)), cbind(samDPcf[i,])), type='l', col=barcols[i], lwd=0.5, lend=3)
      graphics::points(cbind(cbind(samDScf[i,]), cbind(rep(samDP[i],2))), type='l', col=barcols[i], lwd=0.5, lend=3)
    }
    graphics::points(samDS, samDP, pch=pchpt, col=barcols, bg=bgpt, cex=0.8)
  }
  graphics::points(rep(meanDS, 2), DPcf, type='l', col=colmean, lwd=2, lend=3)
  graphics::points(DScf, rep(meanDP, 2), type='l', col=colmean, lwd=2, lend=3)
  if (addlegend & length(x$gp) > 0) graphics::legend('topleft', pch=pchpt, col=colci.gp, levels(x$gp))
  graphics::box()
}
