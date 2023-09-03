#' A Bivariate Plot of Live-Dead Offsets in Alpha Diversity and Evenness
#'
#' AlphaPlot function generates a bivariate plot of offsets in alpha diversity (DELTA S; x-axis)
#' and evenness (DELTA PIE; y-axis) for within-site pairwise comparisons of live and dead samples.
#'
#' @details AlphaPlot function uses the object of the class FidelityDiv
#' by \code{\link{FidelityDiv}} function to produce a cross-plot visualizing
#' live-dead differences in alpha diversity (DELTA S) and evenness (DELTA PIE),
#' including confidence bars for each site. If a grouping factor is provided,
#' symbols and bars are color-coded by levels and group means are plotted
#' in addition to grand mean for all sites. Other graphic arguments
#' applicable to function \code{\link[graphics]{plot.default}} can be specified.
#'
#' @param x An object of the class 'FidelityDiv' returned by  \code{\link{FidelityDiv}} function
#'
#' @param pch An integer defining symbol type (default=21). Filled symbols (21-25) should be used.
#'
#' @param col A character string or number defining color of symbols for individual sites
#'  (default = 'black'). This argument applies only if a grouping factor 'gp' is not provided.
#'
#' @param bgpt A character string or number defining background color of symbols for individual
#'  sites (default = 'white').
#'
#' @param transp Numerical (default=1): defines transparency of symbols and associated
#'  confidence intervals. Allowed values range from 0 (transparent) to 1 (opaque).
#'  Note: transp=0 will suppress plotting of symbols and bars.
#'
#' @param cex A numerical value (default=0.8) defining expansion factor for symbols.
#'
#' @param cf.bars Logical (default=TRUE) to indicate if confidence bars should be plotted.
#'
#' @param cf.col A character string or a number defining color of confidence bars for individual
#'  sites (default = 'black'). This argument applies only if a grouping factor 'gp' is not provided.
#'
#' @param cf.lwd A numerical value defining line width for confidence bars
#' (default = 1).
#'
#' @param colmean A character string or number defining color of bars representing
#'  mean estimates for all sites (default = 'black')
#'
#' @param col.gp A character string or number defining colors for group (default = 1:#levels)
#'  for symbols and bars for sites (applicable when 'gp' factor is provided). The number of colors
#'  must match number of levels in a grouping factor
#'
#' @param mean.lwd A numerical value defining line width for confidence bars for group means
#' (default = 2).
#'
#' @param transp_mean Numerical (default=1): defines transparency of means and associated
#' confidence intervals. Allowed values range from 0 (transparent = invisible)
#' to 1 (opaque). Note: transp=1 will suppress plotting of symbols and bars for group means.
#'
#' @param xlab A character string defining x axis label (default=bquote(Delta['S']).
#'
#' @param ylab A character string defining y axis label (default=bquote(Delta['PIE']).
#'
#' @param axes Logical (default = TRUE): plots or suppresses axes
#'
#' @param xmax A positive numerical value defining x axis range (-xmax, xmax)
#'
#' @param ymax A positive numerical value defining y axis range (-ymax, ymax)
#'
#' @param legend.cex A numerical value (default = 1) defining text size for legend.
#'
#' @param addlegend Logical (default=TRUE): adds legend to the plot, if 'gp' factor is provided
#'
#' @param ... Other graphic arguments applicable to function \code{\link[graphics]{plot.default}}
#'
#' @return A bivariate plot.
#'
#'
#' @examples
#'
#' my.fid1 <- FidelityDiv(FidData$live, FidData$dead, iter=100, n.filters=50)
#' AlphaPlot(my.fid1)
#' my.fid2 <- FidelityDiv(FidData$live[6:9,], FidData$dead[6:9,], FidData$habitat[6:9],
#'                        n.filters=20, iter=100, CI=0.95)
#' AlphaPlot(my.fid2, col.gp=c('forestgreen', 'coral1'), bgpt='beige', mean.lwd=5, transp_mean=0.5)
#'
#' @export

AlphaPlot <- function(x, pch = 21, col = 'black', bgpt = 'white', transp = 1, cex = 0.8,
                      cf.bars = TRUE, cf.col = 'black', cf.lwd = 1, colmean='black',
                      col.gp = NULL, mean.lwd = 2, transp_mean = 1, xlab = bquote(Delta[' S']),
                      ylab = bquote(Delta[' PIE']), axes = TRUE, xmax = NULL, ymax = NULL,
                      legend.cex=1, addlegend=TRUE, ...) {

  if (!('FidelityDiv' %in% class(x))) stop('object of the class "FidelityDiv" is required')
  if (length(xmax) == 0) xmax <- 1.05 * max(c(abs(x$x[,3]), abs(x$x[,4])))
  if (length(ymax) == 0) ymax <- 1.05 * max(c(abs(x$x[,3]), abs(x$x[,4])))
  graphics::plot(x$x[,2], x$y[,2], xlim=c(-xmax, xmax), ylim=c(-ymax, ymax), type='n', las=1,
       xlab=xlab, ylab=ylab, cex.lab=1.5, axes = axes, ...)
  graphics::abline(h=0, v=0, lwd=1, col='gray')
  if (length(x$gp)==0) {
    if (cf.bars) {
     for(i in 1:length(x$x[,2])) {
     graphics::points(cbind(cbind(rep(x$x[i,2], 2)), cbind(x$y[i,3:4])), type='l',
                      col=grDevices::adjustcolor(cf.col, transp), lwd=cf.lwd, lend=3)
     graphics::points(cbind(cbind(x$x[i,3:4]), cbind(rep(x$y[i,2],2))), type='l',
                      col=grDevices::adjustcolor(cf.col, transp), lwd=cf.lwd, lend=3)
     }
    }
    graphics::points(x$x[,2], x$y[,2], pch=pch, col=grDevices::adjustcolor(col, transp),
                     bg=bgpt, cex=cex)
  }
  if (length(x$gp)>0) {
    ifelse(length(col.gp) > 0, col.gp <- col.gp, col.gp <- 1:(length(levels(x$gp))))
    barcols <- col.gp[x$gp]
    if (cf.bars) {
     for(i in 1:length(x$x[,2])) {
      graphics::points(cbind(cbind(rep(x$x[i,2], 2)), cbind(x$y[i,3:4])), type='l',
                       col=grDevices::adjustcolor(barcols[i], transp), lwd=cf.lwd, lend=3)
      graphics::points(cbind(cbind(x$x[i,3:4]), cbind(rep(x$y[i,2],2))), type='l',
                       col=grDevices::adjustcolor(barcols[i], transp), lwd=cf.lwd, lend=3)
     }
    }
    graphics::points(x$x[,2], x$y[,2], pch=pch,
                     col=grDevices::adjustcolor(barcols, transp), bg=bgpt, cex=cex)
  }

if (nrow(x$x) > 1) {

  if (length(x$xgp) == 0) {
   graphics::points(rep(x$xmean[1], 2), x$ymean[2:3], type='l',
                    col=grDevices::adjustcolor(colmean, transp_mean), lwd=mean.lwd, lend=3)
   graphics::points(x$xmean[2:3], rep(x$ymean[1], 2), type='l',
                    col=grDevices::adjustcolor(colmean, transp_mean), lwd=mean.lwd, lend=3)
  }

  if (length(x$xgp)>0) {
   for(i in 1:nrow(x$xgp)) {
    graphics::points(rep(x$xgp[i,2], 2), x$ygp[i,3:4], type='l',
                    col=grDevices::adjustcolor(col.gp[i], transp_mean), lwd=mean.lwd, lend=3)
    graphics::points(x$xgp[i,3:4], rep(x$ygp[i,2], 2), type='l',
                    col=grDevices::adjustcolor(col.gp[i], transp_mean), lwd=mean.lwd, lend=3)
   }
  }

}

  if (addlegend & length(x$gp) > 0) graphics::legend('topleft', pch=pch, pt.bg=bgpt, cex=legend.cex,
                                                     col=col.gp, levels(x$gp), pt.cex=legend.cex)
  graphics::box()
}
