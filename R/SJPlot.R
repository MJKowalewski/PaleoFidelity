#' A bivariate fidelity plot
#'
#' SJPlot function generates a bivariate plot of a correlation measure (x axis) versus
#' a similarity measure (y axis) for all live-dead pairwise comparisons.
#'
#' @details SJPlot function is designed to use the output of FidelityEst function to produce
#' a classic fidelity plot (see Kidwell 2007). Specifically, if default arguments
#' for fidelity measures ('Spearman' and 'Chao') are used in FidelityEst function,
#' a Spearman vs. Jaccard-Chao fidelity plot (as in Kidwell, 2007) is produced.
#' If a grouping factor is provided, symbols are color-coded by levels and group means are plotted.
#' Bivariate distributions produced by resampling models can also be included.
#'
#' NOTE: This is a simple wrap-up of plot function, including some of its common graphic arguments.
#' It allows for quick exploratory plots and should be readily editable to derive more customized plots.
#'
#' @param x An object (a list) returned by FidelityEst function.
#'
#' @param bubble Logical (default=TRUE): to produce a bubble plot with symbols scaled by N min
#' (the number of obeservations in the smaller of the two (live vs. dead) compared samples).
#'
#' @param xlim A vector with two numerical values representing x axis limits (default = c(-1, 1)).
#'
#' @param ylim A vector with two numerical values representing y axis limits (default = c(0, 1)).
#'
#' @param trans A numeric value (default = 0.3) defining transparency of background fill for symbols.
#'
#' @param cex A numeric value (default = 1) defining symbol size (applicable if Bubble = FALSE).
#'
#' @param pch An integer or a single character (default = 21) specifying symbol type.
#'
#' @param col Color name (default = 'black') defining symbol color, applicable when 'gp' factor
#'  was provided in FidelityEst function.
#'
#' @param gpcol Color names (default = 1:#levels) defining colors for sample groups, applicable
#'  when 'gp' factor was provided in FidelityEst function. If custom colors are provided, the number
#'  of colors must match number of levels in a 'gp' factor.
#'
#' @param pch2 An integer or a single character (default = 21) specifying symbol type
#' for grand or group means.
#'
#' @param PF Logical (deafault=FALSE): to add scatterplot of resampled fidelity estimates
#'  for the null model postulating perfect fidelity.
#'
#' @return A single bivariate plot produce by plot function.
#'
#'
#' @examples
#'
#' out1 <- FidelityEst(FidData$live, FidData$dead, n.filters=30, t.filters=1)
#' SJPlot(out1)
#'
#' @export

SJPlot <- function(x, bubble=TRUE, xlim=c(-1, 1), ylim=c(0, 1), trans=0.3, cex=1,
                        pch=21, col='black', gpcol=NULL, pch2='+', PF=FALSE)
{
  graphics::plot(x$x, x$y, type='n', xlim=xlim, ylim=ylim,
                 xlab=x$measures[1], ylab=x$measures[2], las=1)

  graphics::abline(h = 0.5, v = 0, lwd = 1.5, col = 'darkgray')

  if (PF) {
    graphics::points(x$PF.dist[,1], x$PF.dist[,2], pch='.', col='coral3')
    graphics::points(mean(x$PF.dist[,1]), mean(x$PF.dist[,2]), pch=21, col='coral3', bg='white', cex=1)
  }

  if (bubble) {
    cexR <- apply(cbind(rowSums(x$live), rowSums(x$dead)), 1, min)
    if (max(cexR) - min(cexR) == 0) cex=cex
    else cex <- 2.5*(0.3 + (cexR - min(cexR)) / (max(cexR) - min(cexR)))
    graphics::legend('topleft', pch=pch, col=col, pt.cex=c(max(cex), min(cex)),
                     as.character(c(max(cexR),min(cexR))), title = 'N min')
  }
  if (length(x$gp) > 0) {
    ifelse(length(gpcol) == 0, gpcol <- 1:length(levels(x$gp)), gpcol <- gpcol)
    graphics::points(x$x, x$y, pch=pch, bg=grDevices::adjustcolor(gpcol, trans)[x$gp],
                     col=gpcol[x$gp], cex=cex)
    graphics::points(x$observed.means[-1,1:2], pch=pch2, col=gpcol, cex=2)
    graphics::legend('bottomleft', pch=pch, col=gpcol,
                     pt.bg=grDevices::adjustcolor(gpcol, trans), pt.cex=1,
                     levels(x$gp), title = 'groups')
  }
  else {
    graphics::points(x[[1]], x[[2]], pch=pch, bg=grDevices::adjustcolor(col, trans), col=col, cex=cex)
    graphics::points(rbind(x$observed.means[1:2]), pch=pch2, col=col, cex=2)
  }
}

