#' A Fidelity Plot
#'
#' SJPlot function produces a bivariate plot of a correlation measure (x axis)
#' versus similarity measure (y axis) to depict live-dead fidelity at one or more
#' sites.
#'
#' @details SJPlot function uses the object of the class 'FidelityEst', produced by
#' \code{\link{FidelityEst}} function, to generate a live-dead fidelity plot.
#' Specifically, if default arguments for fidelity measures are used
#' in \code{\link{FidelityEst}} function, a Spearman vs. Jaccard-Chao fidelity plot
#' (as in Kidwell, 2007) is produced. If a grouping factor is provided, symbols are
#' color-coded by levels and group means are plotted. Bivariate distributions
#' produced by resampling models can also be included.
#'
#' NOTE: This function utilizes \code{\link[graphics]{plot.default}} function,
#' including some of its common graphic arguments. It allows to explore visually
#' various alternative estimates of fidelity and is readily editable to derive
#' more customized plots.
#'
#' @param x An object of the class 'FidelityEst' returned by \code{\link{FidelityEst}} function.
#'
#' @param bubble Logical (default=TRUE): to produce a bubble plot with symbols scaled by N-min
#' (the number of observations in the smaller of the two (live vs. dead) compared samples).
#'
#' @param xlim A vector with two numerical values defining x axis limits (default = c(-1, 1)).
#'
#' @param ylim A vector with two numerical values defining y axis limits (default = c(0, 1)).
#'
#' @param trans A numerical value (default = 0.3) defining transparency of background fill for symbols.
#'
#' @param cex A numerical value (default = 1) defining symbol size (applicable if Bubble = FALSE).
#'
#' @param cex.mean A numerical value (default = 2) defining symbol size for means
#' across all samples or sample groups.
#'
#' @param legend.cex A numerical value (default = 1) defining text size for legend.
#'
#' @param axes Logical (default = TRUE): to determine if axes should be plotted
#'
#' @param pch An integer or a single character (default = 21) specifying symbol type.
#'
#' @param col A character string (default = 'black') defining symbol color, applicable when 'gp' factor
#'  was provided in FidelityEst function.
#'
#' @param gpcol Color names (default=1:length(levels(x$gp))) defining colors for sample groups,
#' applicable when 'gp' factor was provided in \code{\link{FidelityEst}} function.
#' If custom colors are provided, the number of colors must match number of levels
#' in a 'gp' factor.
#'
#' @param pch2 An integer or a single character (default=21) specifying symbol type
#' for grand mean or group means.
#'
#' @param PF Logical (default=FALSE): to add a scatter plot of resampled fidelity
#'  estimates for the null model postulating perfect fidelity.
#'
#' @param CI Logical (default=TRUE): adds confidence intervals to corrected fidelity estimates
#'
#' @param adjF Logical (default=TRUE): plots adjusted fidelity estimates
#'
#' @param ssF Logical (default=FALSE): plots sample-standardized fidelity estimates
#'
#' @param unadjF Logical (default=FALSE): plots unadjusted fidelity estimates
#'
#' @param addlegend Logical (default=TRUE): adds legend to the plot, if 'gp' factor is provided.
#'
#' @param addinfo Logical (default=TRUE): prints parameter values above the plot
#'
#' @param addbubble Logical (default=TRUE): prints legend relating bubble sizes to
#' N-min, which denotes the number of specimens in the smaller of the two compared
#' samples (usually, 'live' sample is the smaller of the two)
#'
#' @param info.y Numeric (default = 0.5): specifies how far above the plot should
#' "addinfo" legend (see above) be plotted. This is "line" parameter value
#' in \code{\link[graphics]{mtext}} function)
#'
#' @param bubble.y Numeric (default = - 0.15): specifies how far above the plot should
#' "addbubble" legend (see above) be plotted. This is "inset" parameter value
#' in \code{\link[graphics]{legend}} function)
#'
#' @param xlab Character string that provides x-axis label. The default value is inherited
#' from FidelityEst object based on a correlation measure used there.
#'
#' @param ylab Character string that provides y-axis label. The default value is inherited
#' from FidelityEst object based on a similarity measure used there.
#'
#' @return A single bivariate plot.
#'
#'
#' @examples
#'
#' out1 <- FidelityEst(FidData$live[6:9,], FidData$dead[6:9,], FidData$habitat[6:9],
#'                    n.filters=30, iter=100, t.filters=1)
#' SJPlot(out1)
#'
#' @export
#'
#' @references Kidwell, S.M., 2007, Discordance between living and death assemblages
#' as evidence for anthropogenic ecological change. Proc Natl Acad Sci USA 104(45): 17701–17706.


SJPlot <- function(x, bubble = TRUE, xlim = c(-1, 1), ylim = c(0, 1), trans = 0.3,
                   cex = 1, cex.mean = 2, legend.cex = 0.8, axes = T, pch = 21,
                   col = 'black', gpcol = NULL, pch2 = '+', PF = TRUE, CI = TRUE,
                   adjF = TRUE, ssF = FALSE, unadjF = FALSE, addlegend = TRUE,
                   addinfo = TRUE, addbubble = TRUE, info.y = 0.5, bubble.y = -0.15,
                   xlab = NULL, ylab = NULL)
{

  if (!('FidelityEst' %in% class(x))) stop('object of the class "FidelityEst" is required')
  if (length(xlab) == 0) xlab <- x$values$measures[1]
  if (length(ylab) == 0) ylab <- x$values$measures[2]

  graphics::plot(x$xc[,1], x$yc[,1], type='n', xlim=xlim, ylim=ylim, las=1,
                 xlab=xlab, ylab=ylab, axes=axes)

  graphics::abline(h = 0.5, v = 0, lwd = 1.5, col = 'darkgray')

  if (PF) {
    if (length(unlist(x$x)) > 1) {
    for (i in 1:ncol(x$x.pf.dist)) {
     graphics::points(x$x.pf.dist[,i], x$y.pf.dist[,i], pch='.', col='coral3')
     graphics::points(mean(x$x.pf.dist[,i]), mean(x$y.pf.dist[,i]), pch=21, col='coral3',
                     bg='white', cex=1)
     }
    }
    else {
     graphics::points(x$x.pf.dist, x$y.pf.dist, pch='.', col='coral3')
     graphics::points(mean(x$x.pf.dist), mean(x$y.pf.dist), pch=21, col='coral3',
                       bg='white', cex=1)
    }
   }

  if (bubble & length(unlist(x$x)) > 1) {
    cexR <- apply(cbind(rowSums(x$live), rowSums(x$dead)), 1, min)
    if (max(cexR) - min(cexR) == 0) cex=cex
    else cex <- 2 * (0.3 + (cexR - min(cexR)) / (max(cexR) - min(cexR)))
    if (addbubble) {
    graphics::legend('topleft', pch = pch, col = col, pt.cex = c(max(cex), min(cex)),
                     cex = legend.cex, as.character(c(max(cexR), min(cexR))),
                     title = '', xpd = NA, inset = c(0, bubble.y), ncol = 2, bty = 'n')
    }
  }
  if (length(x$gp) > 0 & length(unlist(x$x)) > 1) {
    ifelse(length(gpcol) == 0, gpcol <- 1:length(levels(x$gp)), gpcol <- gpcol)
    if (CI & adjF) graphics::arrows(x$xc[,2], x$yc[,1], x$xc[,3], x$yc[,1],
                                    length=0, col=gpcol[x$gp], lwd=0.3)
    if (CI & adjF) graphics::arrows(x$xc[,1], x$yc[,2], x$xc[,1], x$yc[,3],
                                    length=0, col=gpcol[x$gp], lwd=0.3)
    if (CI & ssF) graphics::arrows(x$xs[,2], x$ys[,1], x$xs[,3], x$ys[,1],
                                    length=0, col=gpcol[x$gp], lwd=0.3)
    if (CI & ssF) graphics::arrows(x$xs[,1], x$ys[,2], x$xs[,1], x$ys[,3],
                                    length=0, col=gpcol[x$gp], lwd=0.3)
    if (adjF) graphics::points(x$xc[,1], x$yc[,1], pch=pch, col=gpcol[x$gp], cex=cex,
                               bg=grDevices::adjustcolor(gpcol, trans)[x$gp])
    if (adjF & length(unlist(x$x)) > 1) graphics::points(x$xc.stats[-1,1], x$yc.stats[-1,1],
                                        pch=21, col=gpcol, bg='white', cex=cex.mean)
    if (ssF) graphics::points(x$xs[,1], x$ys[,1], pch=23, col=gpcol[x$gp], cex=cex,
                                 bg=grDevices::adjustcolor(gpcol, trans)[x$gp])
    if (ssF) graphics::points(x$xs.stats[-1,1], x$ys.stats[-1,1],
                                        pch=21, col=gpcol, bg='white', cex=cex.mean)
    if (unadjF) graphics::points(unlist(x$x), unlist(x$y), pch=pch, col=gpcol[x$gp], cex=cex,
                                 bg=grDevices::adjustcolor(gpcol, trans)[x$gp])
    if (unadjF) graphics::points(x$x.stats[-1,1], x$y.stats[-1,1], pch=pch,
                              col=gpcol, bg='white', cex=cex.mean)
    if (addlegend) {
      graphics::legend('topleft', pch=pch, col=gpcol, cex=legend.cex,
                     pt.bg=grDevices::adjustcolor(gpcol, trans), pt.cex=legend.cex,
                     levels(x$gp), title = 'groups')
     }
   }
  else {
    if (adjF) graphics::arrows(x$xc[,1], x$yc[,2], x$xc[,1], x$yc[,3],
                               length=0, col=col, lwd=0.1)
    if (adjF) graphics::arrows(x$xc[,2], x$yc[,1], x$xc[,3], x$yc[,1],
                               length=0, col=col, lwd=0.1)
    if (adjF) graphics::points(x$xc[,1], x$yc[,1], pch=pch, col=col, cex=cex,
                               bg=grDevices::adjustcolor(col, trans))
    if (adjF & length(unlist(x$x)) > 1) graphics::points(x$xc.stats[1,1], x$yc.stats[1,1],
                               pch=pch2, col=col, bg='white', cex=cex.mean)
    if (ssF) graphics::arrows(x$xs[,1], x$ys[,2], x$xs[,1], x$ys[,3],
                               length=0, col=col, lwd=0.1)
    if (ssF) graphics::arrows(x$xs[,2], x$ys[,1], x$xs[,3], x$ys[,1],
                               length=0, col=col, lwd=0.1)
    if (ssF) graphics::points(x$xs[,1], x$ys[,1], pch=pch, col=col, cex=cex,
                              bg=grDevices::adjustcolor(col, trans))
    if (ssF & length(unlist(x$x)) > 1) graphics::points(x$xs.stats[1,1], x$ys.stats[1,1],
                               pch=pch2, col=col, bg='white', cex=cex.mean)
    if (unadjF) graphics::points(unlist(x$x), unlist(x$y), pch=pch, col=col, cex=cex,
                                 bg=grDevices::adjustcolor(col, trans))
    if (unadjF) graphics::points(mean(unlist(x$x)), mean(unlist(x$y)),
                               pch=pch, col=col, bg='white', cex=cex.mean)

  }
 if (addinfo) {
if (adjF & !ssF)   graphics::mtext(side = 3, line = info.y, cex = legend.cex,
                   paste('transform=', x$values$data.transf,
                               '   ', 'PF model iter=', x$values$PFiter, sep=''))
if (ssF & !adjF)   graphics::mtext(side = 3, line = info.y, cex = legend.cex,
                              paste('transform=', x$values$data.transf,
                                    '   ', 'subsampling iter=', x$values$SSIter, sep=''))
if (ssF & adjF)   graphics::mtext(side = 3, line = info.y, cex = legend.cex,
                                      paste('transform=', x$values$data.transf,
                                            '   ','PF model iter=', x$values$PFiter,'   ',
                                            'subsampling iter=', x$values$SSIter, sep=''))

 }
}
