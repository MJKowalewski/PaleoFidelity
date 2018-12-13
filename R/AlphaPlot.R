#' AlphaPlot (a cross plot of live-dead offsets in alpha diversity and evenness)
#'
#' AlphaPlot function generates a bivariate plot of an offset in alpha diversity (DELTA S)
#' on x axis versus an offset in evenness (DELTA PIE) for all live-dead pairwise comparisons.
#'
#' @details AlphaPlot function uses the output of FidelityDiv function to produce
#' a cross plot visualizing differences in estimates of alpha diversty and evenness between
#' pairs of sympatric live and dead samples. In its default form, an approch proposed by
#' Olszewski and Kidwell (2007) is used. The offset in alpha diversity is measured as the difference
#' between natural logarithms of sample standardized richness log(R-dead) - log(R-live). The offset
#' in evenness is measured as the difference between Hurlbert's PIE (PIE-dead - PIE-live).
#' Confidence bars depitcting user-specified confidence intervals are also plotted.
#' If a grouping factor is provided, symbols and bars are color-coded by levels and
#' group means are plotted in addition to grand mean for all sites.
#'
#'
#' @param x An object (a list) returned by FidelityDiv function
#'
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
#' @param cf.bars Logical (default=TRUE) to indicate if confidence bars should be plotted
#'
#' @param cex a numerical value (default=0.8) defining expansion factor for symbols
#'
#' @param pchpt A number defining symbol type (default = 21) (filled symbols (21-25) should be used)
#'
#' @param colmean A character string or number defining color of bars representing
#' mean estimates for all sites (default = 'white')
#'
#' @param col.gp A character string or number defining colors for group (default = 1:#levels)
#' for symbols and bars for sites (applicable when 'gp' factor is provided). The number of colors
#' must match number of levels in a grouping factor
#'
#' @param addlegend Logical (default=TRUE): add legend
#' Note: The legend is plotted only when 'gp' factor provided
#'
#' @return A single bivariate plot produce by plot function
#'
#'
#' @examples
#'
#' my.fid1 <- FidelityDiv(FidData$live, FidData$dead, n.filters=30)
#' AlphaPlot(my.fid1)
#' my.fid2 <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, n.filters=50, iter=1000, CI=0.95)
#' AlphaPlot(my.fid2, col.gp=c('forestgreen', 'coral1'), bgpt='beige')
#'
#' @export

AlphaPlot <- function(x, colpt='gray', colci='gray', bgpt='white', cf.bars=TRUE, cex=0.8,
                      pchpt=21, colmean='black', col.gp=NULL, addlegend=TRUE) {
  xmax <- 1.2 * max(abs(x$x[,2]))
  ymax <- 1.2 * max(abs(x$y[,2]))
  graphics::plot(x$x[,2], x$y[,2], xlim=c(-xmax, xmax), ylim=c(-ymax, ymax), type='n', las=1,
       xlab=bquote(Delta['S']), ylab=bquote(Delta['PIE']), cex.lab=1.5)
  graphics::abline(h=0, v=0, lwd=1, col='gray')
  if (length(x$gp)==0) {
    if (cf.bars) {
     for(i in 1:length(x$x[,2])) {
     graphics::points(cbind(cbind(rep(x$x[i,2], 2)), cbind(x$y[i,3:4])), type='l', col=colci, lwd=0.5, lend=3)
     graphics::points(cbind(cbind(x$x[i,3:4]), cbind(rep(x$y[i,2],2))), type='l', col=colci, lwd=0.5, lend=3)
     }
    }
    graphics::points(x$x[,2], x$y[,2], pch=pchpt, col=colpt, bg=bgpt, cex=cex)
  }
  if (length(x$gp)>0) {
    ifelse(length(col.gp) > 0, col.gp <- col.gp, col.gp <- 1:(length(levels(x$gp))))
    barcols <- col.gp[x$gp]
    if (cf.bars) {
     for(i in 1:length(x$x[,2])) {
      graphics::points(cbind(cbind(rep(x$x[i,2], 2)), cbind(x$y[i,3:4])), type='l', col=barcols[i], lwd=0.5, lend=3)
      graphics::points(cbind(cbind(x$x[i,3:4]), cbind(rep(x$y[i,2],2))), type='l', col=barcols[i], lwd=0.5, lend=3)
     }
    }
    graphics::points(x$x[,2], x$y[,2], pch=pchpt, col=barcols, bg=bgpt, cex=cex)
  }

  if (length(x$xgp) == 0) {
   graphics::points(rep(x$xmean[1], 2), x$ymean[2:3], type='l', col=colmean, lwd=4, lend=3)
   graphics::points(x$xmean[2:3], rep(x$ymean[1], 2), type='l', col=colmean, lwd=4, lend=3)
  }

  if (length(x$xgp)>0) {
   for(i in 1:nrow(x$xgp)) {
   graphics::points(rep(x$xgp[i,2], 2), x$ygp[i,3:4], type='l', col=col.gp[i], lwd=4, lend=3)
   graphics::points(x$xgp[i,3:4], rep(x$ygp[i,2], 2), type='l', col=col.gp[i], lwd=4, lend=3)
   }
  }

  if (addlegend & length(x$gp) > 0) graphics::legend('topleft', pch=pchpt, pt.bg=bgpt,
                                                     col=col.gp, levels(x$gp), pt.cex=cex)
  graphics::box()
}
