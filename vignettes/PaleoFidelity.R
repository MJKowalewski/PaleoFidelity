## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installing PaleoFidelity package, eval = FALSE---------------------------
#  install.packages('devtools')
#  library(devtools)
#  devtools::install_github('mjkowalewski/PaleoFidelity', build_vignettes = TRUE)
#  library(PaleoFidelity)

## ----data example-------------------------------------------------------------
library(PaleoFidelity)
str(FidData) # check the structure of the dataset

## ----fidelity summary function------------------------------------------------
FidelitySummary(live = FidData$live, dead = FidData$dead, gp = FidData$habitat, report = TRUE)

## ----fidelity summary function part 2-----------------------------------------
FidelitySummary(live = FidData$live, dead = FidData$dead, gp = FidData$habitat,
                report = TRUE, n.filters = 30)

## ----fidelity summary function part 3-----------------------------------------
FidelitySummary(live = FidData$live, dead = FidData$dead, gp = FidData$habitat, report = TRUE, n.filters = 100)

## ----live-dead plot, fig.width=7, fig.height=6--------------------------------
par(mar=c(3, 7, 0.5, 7))
rep1 <- LDPlot(live = colSums(FidData$live),
       dead = colSums(FidData$dead),
       tax.names = colnames(FidData$live), toplimit = 20,
       cor.measure = 'spearman', report = TRUE, iter = 1000)

## ----LD comparison------------------------------------------------------------
rep1[1:5]

## ----live-dead model, fig.width=7, fig.height=3.5-----------------------------
par(mar=c(4, 4, 0.5, 0.5))
hist(rep1$randomized.r[,2], breaks=seq(-1,1,0.05), main='',
     las=1, xlab=bquote('Spearman' ~ italic(rho)))
arrows(rep1$cor.coeff[2], 100, rep1$cor.coeff[2], 10,
       length=0.1, lwd=2)

## ----fidelity estimates-------------------------------------------------------
out1 <- FidelityEst(live = FidData$live, dead = FidData$dead,
                    gp = FidData$habitat,
                    n.filters = 30, iter = 499)
str(out1)

## ----fidelity estimates outputs-----------------------------------------------
out1$xc # adjusted correlation measure summary
out1$yc # adjusted similarity measure summary

## ----fidelity estimates outputs: sample-standardized--------------------------
out1$xs # sample-standardized correlation measure summary
out1$ys # sample-standardized similarity measure summary

## ----classic fidelity plot, fig.width=7, fig.height=4-------------------------
par(mar = c(4, 4, 0.5, 0.5))
SJPlot(out1, gpcol = c('aquamarine3', 'coral3'), cex.legend = 0.8)

## ----classic fidelity plot 2, fig.width=7.5, fig.height=4---------------------
par(mar = c(4, 4, 0.5, 0.5))
SJPlot(out1, gpcol = c('aquamarine3', 'coral3'), bubble = F, unadj = T, adjF = F, cex.legend = 0.8)

## ----alpha diversity----------------------------------------------------------
out3 <- FidelityDiv(FidData$live, FidData$dead, iter=1000)
out3$x
out3$y

## ----alpha diversity 2--------------------------------------------------------
out4 <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, iter=1000)
out4$xmean
out4$ymean
out4$xgp
out4$ygp
out4$p.values
out4$p.gps

## ----plot alpha 2, fig.width=7, fig.height=4----------------------------------
out3 <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, CI = 0.95, iter = 1000)
par(mar = c(4, 4.5, 0.5, 0.5))
AlphaPlot(out3, col.gp = c('aquamarine3', 'coral3'), bgpt = 'beige', pch = 22, legend.cex = 0.8)

