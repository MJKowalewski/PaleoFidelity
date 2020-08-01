#' A Comparative Live-Dead Barplot
#'
#' LDPlot function generates a comparative live-dead plot to visualize
#' congruence or discordance in relative abundance and rank order of the most
#' common species, genera, functional groups or other enumerated variables.
#'
#' @details LDPlot function produces a plot that compares barplots for live and
#' dead datasets for the top "n" most common species/variables (e.g.,
#' Kowalewski et al. 2003). These L-D comparisons can apply to single sites
#' or data pooled across multiple sites. Because names of species/variables
#' can vary in length and because the number of plotted variables can span
#' a wide range of values, margin widths for plots and cex parameter
#' may need to be customized on a case by case basis.
#'
#' @param live A vector of integers with counts of live specimens by species (or other units)
#'
#' @param dead A vector of integers with counts of dead specimens by species (or other units)
#'
#' @param tax.names A vector with a list of names of species (or other units used as variables)
#'
#' @param toplimit A numerical value (default = 10) defining the number of top
#' species (or other units) to be plotted
#'
#' @param barwidth A numerical value (default = 250 / toplimit) defining bar width
#' (more precisely, the thickness of \code{\link[graphics]{lines}} used to
#' represent bars)
#'
#' @param col1 A character string (default = 'black') defining the color of bars
#' for species/variables shared by live and dead data
#'
#' @param col2 A character string (default = 'gray') defining the color of bars
#' for species/variables unique to either live or dead data
#'
#' @param arr.col A character string (default = 'black') defining the color of arrows
#'
#' @param arr.lty An integer or character string specifying type of line used for
#' arrows (default = 1) (this value is passed on to the graphical parameter 'lty')
#'
#' @param arr.lwd A numerical value  specifying the width of arrow lines (default = 1)
#' (this value is passed on to the graphical parameter 'lwd')
#'
#' @param cex.axis A numerical value (default = 0.7).The magnification to be used for
#' x axis annotation relative to the current setting of cex.
#'
#' @param tck A numerical value (default  = -0.02) defining length of tick marks
#' for x axis
#'
#' @param cex.lab A numerical value (default = 0.8) defining font size for
#' x-axis label
#'
#' @param cex.names A numerical value (default = 0.95) defining font size
#' for names of species/variables corresponding to individual bars on the barplot
#'
#' @param cex.label A numerical value (defult = 1) defining font size for 'Live'
#' and 'Dead' titles placed above the barplots
#'
#' @param cex.stat A numerical value (default = 0.9) defining font size
#' for correlation coefficient estimate
#'
#' @param font.names A numerical value (default = 3 = italic) defining font style
#'
#' @param cor.measure A character string (default='spearman') defining correlation measure
#'  (passed on to \code{\link[stats]{cor}} function) used to estimate
#'  live-dead correlations (possible values are 'pearson', 'spearman', 'kendall', 'all').
#'  When 'all' is selected, all three correlation measures are printed on the plot.
#'  Note that the computed statistics include all taxa regardless of how many taxa are
#'  plotted on the chart as defined by 'toplimit' parameter.
#'
#' @return A single plot.
#'
#' @examples
#'
#' temp.par <- par(mar=c(3,6,1,6))
#' LDPlot(live=colSums(FidData$live), dead=colSums(FidData$dead),
#' tax.names=colnames(FidData$live), toplimit=15, barwidth = 21,
#' col1 = 'green2', col2 = 'red4', arr.col = 'green2', arr.lty = 1)
#' par(temp.par)
#'
#' @export
#'
#' @importFrom stats cor
#'
#' @references Kowalewski, M.,	Carroll, M., Casazza, L., Gupta, N., Hannisdal, B.,
#' Hendy, A., Krause, R.A., Jr., Labarbera, M., Lazo, D.G., Messina, C., Puchalski, S.,
#' Rothfus, T.A., Sälgeback, J., Stempien, J., Terry, R.C., Tomašových, A., (2003),
#' Quantitative fidelity of brachiopod-mollusk assemblages from modern subtidal
#' environments of San Juan Islands, USA. Journal of Taphonomy 1: 43-65.


LDPlot <- function(live, dead, tax.names, toplimit = 10, barwidth = 150 / toplimit,
                   col1 = 'black', col2 = 'gray', arr.col = 'black', arr.lty=1,
                   arr.lwd = 1, cex.axis = 0.7, tck = -0.02, cex.lab = 0.8,
                   cex.names = 0.95, cex.label = 1, cex.stat = 0.9, font.names = 3,
                   cor.measure = 'spearman') {

  if (!is.vector(live)) stop('object "live" must be a vector')
  if (!is.vector(dead)) stop('object "dead" must be a vector')
  if (!is.numeric(live)) stop('object "live" must be numeric')
  if (!is.numeric(dead)) stop('object "dead" must be numeric')
  if (length(live) != length(dead))
    stop('objects "live" and "dead" are of different length')
  if (length(live) != length(tax.names))
    stop('object "tax.names" differs in length from objects "live" and "dead"')
  if (sum((live + dead) > 0) < toplimit)
    stop(paste('toplimit = ', toplimit, 'exceeds the total number of non-zero
               taxa/variables =', sum((live + dead) > 0), 'set toplimit <=',
               sum((live + dead) > 0)))

# info on top n live and dead units
  toplive <- (live[order(live, decreasing = T)] / sum(live)) [1:toplimit]
  toplivenames <- tax.names[order(live, decreasing=T)][1:toplimit]
  topdead <- (dead[order(dead, decreasing = T)] / sum(dead)) [1:toplimit]
  topdeadnames <- tax.names[order(dead, decreasing = T)] [1:toplimit]
  dead.live <- as.numeric(topdeadnames %in% toplivenames)
  live.dead <- as.numeric(toplivenames %in% topdeadnames)

# compute correlation coefficients
  if (cor.measure != 'all') corest <- round(stats::cor(live, dead, method=cor.measure), 3)
  if (cor.measure == 'pearson') cor.label <- bquote(italic(r) == .(corest))
  if (cor.measure == 'spearman') cor.label <- bquote(italic(rho) == .(corest))
  if (cor.measure == 'kendall') cor.label <- bquote(italic(tau) == .(corest))
  if (cor.measure == 'all') {
    corest1 <- round(stats::cor(live, dead, method='pearson'), 3)
    cor.label1 <- bquote(italic(r) == .(corest1))
    corest2 <- round(stats::cor(live, dead, method='spearman'), 3)
    cor.label2 <- bquote(italic(rho) == .(corest2))
    corest3 <- round(stats::cor(live, dead, method='kendall'), 3)
    cor.label3 <- bquote(italic(tau) == .(corest3))
    }

# customize info for x axis
  max.x <- ceiling(5*max(c(toplive, topdead)))/5
  xlim2 <- 2 * max.x * 1.25
  if (max.x <= 0.5) {
    my.x.at <- c(seq(0, max.x, 0.1), seq(xlim2-max.x, xlim2, 0.1))
    my.x.lab <- c(seq(0, max.x, 0.1), seq(max.x, 0, -0.1))
  }
  if (max.x > 0.5) {
    my.x.at <- c(seq(0, max.x, 0.2), seq(xlim2-max.x, xlim2, 0.2))
    my.x.lab <- c(seq(0, max.x, 0.2), seq(max.x, 0, -0.2))
  }

# plot LD chart
  graphics::plot(0, 0, type = 'n', ylim = c(toplimit, 0), xlim = c(0, xlim2),
                 xlab = '', ylab = '', axes = F)
   for(i in 1:toplimit) {
    if (live.dead[i] == 0) graphics::lines(c(0, toplive[i]), c(i , i),
                                           lwd = barwidth, lend = 3, col = col2)
    if (live.dead[i] == 1) graphics::lines(c(0, toplive[i]), c(i, i),
                                           lwd = barwidth, lend = 3, col = col1)
    if (dead.live[i] == 0) graphics::lines(c(xlim2 - topdead[i], xlim2), c(i, i),
                                           lwd = barwidth, lend = 3, col = col2)
    if (dead.live[i] == 1) graphics::lines(c(xlim2 - topdead[i], xlim2), c(i, i),
                                           lwd = barwidth, lend = 3, col = col1)
  }
  graphics::axis(2, labels = toplivenames, at = 1:toplimit, lwd = 0, las = 1,
                 padj = 0.5, hadj = 1, font = font.names, cex.axis = cex.names)
  graphics::axis(4, labels = topdeadnames, at = 1:toplimit, lwd = 0, las = 1,
                 padj = 0.5, hadj = 0, font = font.names, cex.axis = cex.names)
  graphics::axis(1, labels = my.x.lab, padj = -1.5, tck = tck, at = my.x.at,
                 cex.axis = cex.axis)
  graphics::lines(c(0.45 * xlim2, 0.55 * xlim2), c(toplimit, toplimit), lwd = 50,
                  col = 'white', lend = 3, xpd = NA)
  graphics::mtext(side = 1, line = 1.5, 'proportion of specimens', cex=cex.lab)
  graphics::text(0.1, -0.5, 'LIVE', cex = cex.label, xpd = NA)
  graphics::text(xlim2 - 0.1, -0.5, 'DEAD', cex = cex.label, xpd = NA)
  if (cor.measure != 'all') graphics::text(0.5 * xlim2, -0.5, cor.label, cex = cex.stat)
  if (cor.measure == 'all') {
    graphics::text(0.5 * xlim2, -1.5, cor.label1, cex = cex.stat, xpd=NA)
    graphics::text(0.5 * xlim2, -0.5, cor.label2, cex = cex.stat, xpd = NA)
    graphics::text(0.5 * xlim2, 0.5, cor.label3, cex = cex.stat, xpd = NA)
    }
  for(i in 1:toplimit) {
    if (dead.live[i] == 1) {
      k <- which(toplivenames == topdeadnames[i])
      graphics::arrows((xlim2 - topdead[i]) - (0.03 * xlim2), i,
                       (0.03 * xlim2) + toplive[k], k, length = 0.1,
                       lwd = arr.lwd, code = 3, col = arr.col, lty = arr.lty)
    }
  }
}
