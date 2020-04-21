#' A Comparative Live-Dead Barplot
#'
#' LDPlot function generates a comparative live-dead plot of the most frequent taxa
#' (or other units/variables describing samples) to highlight conguences and discordances
#' in relative abundance and rank order of most abundant taxa/variables.
#'
#' @details LDPlot function is designed to produce a plot with two barplots comparing
#' relative abundances of the most common taxa in live and dead datasets. This plot
#' is suitable for single live-dead comparisons only (single samples,
#' sets of samples or datasets).
#'
#' NOTE: This function utilizes graphics::plot function, including some of its common
#' graphic arguments.
#'
#' @param live A vector with integer counts of live specimens
#'
#' @param dead A vector with integer counts of dead specimens
#'
#' @param tax.names A vector with a list of names of taxa (or other units
#' used as variables)
#'
#' @param toplimit A numerical value (default = 10) defining the number of top
#' species (or other variables) to be plotted
#'
#' @param barwidth A numerical value (default = 250 / toplimit) defining bar width
#' (or more precisely the thicknes of graphics::lines used to represent bars)
#'
#' @param col1 A character string (default = 'black') defining the color of bars
#' for taxa (or other variables) shared by live and dead data
#'
#' @param col2 A character string (default = 'gray') defining the color of bars
#' for taxa (or other variables) unique to either live or dead data
#'
#' @param arr.col A character string (default = 'red') defining the color of arrows
#'
#' @param arr.lty An integer or character string specifying type of line used for
#' arrows (default = 1) (this value is passed on to the graphical parameter 'lty')
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
#' for names of taxa/variables corresponding to individual bars on the barplot
#'
#' @param cex.label A numerical value (defult = 1) defining font size for 'Live'
#' and 'Dead' titles placed above the barplots
#'
#' @return A single plot.
#'
#' @examples
#'
#' temp.par <- par(mar=c(3,5,2,5))
#' LDPlot(live=colSums(FidData$live), dead=colSums(FidData$dead),
#'        tax.names=colnames(FidData$live), toplimit=15, col1 = 'green',
#'        col2 = 'red4', arr.col = 'black', arr.lty = 3)
#' par(temp.par)
#'
#' @export

LDPlot <- function(live, dead, tax.names, toplimit = 10, barwidth = 150 / toplimit,
                   col1 = 'black', col2 = 'gray', arr.col = 'red', arr.lty=1,
                   cex.axis = 0.7, tck = -0.02, cex.lab = 0.8, cex.names = 0.95,
                   cex.label = 1) {

  if (!is.vector(live)) stop('object "live" is not a vector')
  if (!is.vector(dead)) stop('object "live" is not a vector')
  if (length(live) != length(dead))
    stop('objects "live" and "dead" are of different length')
  if (length(live) != length(tax.names))
    stop('object "tax.names" differs in length from "live" and "dead" objects')

  toplive <- (live[order(live, decreasing = T)] / sum(live)) [1:toplimit]
  toplivenames <- tax.names[order(live, decreasing=T)][1:toplimit]
  topdead <- (dead[order(dead, decreasing = T)] / sum(dead)) [1:toplimit]
  topdeadnames <- tax.names[order(dead, decreasing = T)] [1:toplimit]
  dead.live <- as.numeric(topdeadnames %in% toplivenames)
  live.dead <- as.numeric(toplivenames %in% topdeadnames)
  max.x <- ceiling(10*max(c(toplive, topdead)))/10
  xlim2 <- 2 * max.x * 1.25
  if (max.x <= 0.5) {
    my.x.at <- c(seq(0, max.x, 0.1), seq(xlim2-max.x, xlim2, 0.1))
    my.x.lab <- c(seq(0, max.x, 0.1), seq(max.x, 0, -0.1))
  }
  if (max.x > 0.5) {
    my.x.at <- c(seq(0, max.x, 0.2), seq(xlim2-max.x, xlim2, 0.2))
    my.x.lab <- c(seq(0, max.x, 0.2), seq(max.x, 0, -0.2))
  }

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
                 padj = 0.5, hadj = 1, font = 3, cex.axis = cex.names)
  graphics::axis(4, labels = topdeadnames, at = 1:toplimit, lwd = 0, las = 1,
                 padj = 0.5, hadj = 0, font = 3, cex.axis = cex.names)
  graphics::axis(1, labels = my.x.lab, padj = -1.5, tck = tck, at = my.x.at,
                 cex.axis = cex.axis)
  graphics::lines(c(0.42 * xlim2, 0.58 * xlim2), c(toplimit, toplimit), lwd = 50,
                  col = 'white', lend = 3, xpd = NA)
  graphics::mtext(side = 1, line = 1.5, 'proportion of specimens', cex=cex.lab)
  graphics::mtext(side = 3, line = 0, adj = 0, 'LIVE', cex=cex.label)
  graphics::mtext(side = 3, line = 0, adj = 1, 'DEAD', cex=cex.label)
  for(i in 1:toplimit) {
    if (dead.live[i] == 1) {
      k <- which(toplivenames == topdeadnames[i])
      graphics::arrows((xlim2 - topdead[i]) - (0.03 * xlim2), i,
                       (0.03 * xlim2) + toplive[k], k, length = 0.1,
                       lwd = 1, code = 3, col = arr.col, lty = arr.lty)
    }
  }
}
