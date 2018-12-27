## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data example--------------------------------------------------------
library(PaleoFidelity)
str(FidData) # check the structure of the example dataset

## ----fidelity summary function-------------------------------------------
FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, report=TRUE)

## ----fidelity summary function part 2------------------------------------
FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat,
                report=TRUE, n.filters=30)

## ----fidelity summary function part 3------------------------------------
FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat,
                report=TRUE, n.filters=100)

## ----fidelity estimates--------------------------------------------------
out1 <- FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, n.filters=30, dbzero=F, iter=99)
out1$observed.means

## ----fidelity estimates part 2-------------------------------------------
out1$PF.stats

## ----fidelity estimates part 3-------------------------------------------
out1$gp.prob

## ----classic fidelity plot, fig.width=7.5, fig.height=4------------------
par(mar=c(4, 4, 0.5, 0.5))
SJPlot(out1, gpcol=c('aquamarine3', 'coral3'), legend.cex=0.8)

## ----classic fidelity plot 2, fig.width=7.5, fig.height=4----------------
par(mar=c(4, 4, 0.5, 0.5))
SJPlot(out1, gpcol=c('aquamarine3', 'coral3'), bubble=F, legend.cex=0.8)

## ----alpha 1-------------------------------------------------------------
out2 <- FidelityDiv(FidData$live, FidData$dead, iter=1000)
str(out2)

## ----plot alpha 1, fig.width=7.5, fig.height=4---------------------------
par(mar=c(4, 4.5, 0.5, 0.5))
AlphaPlot(out2, legend.cex=0.8)

## ----alpha 2-------------------------------------------------------------
out3 <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, CI=0.95, iter=1000)
str(out3)

## ----plot alpha 2, fig.width=7.5, fig.height=4---------------------------
par(mar=c(4, 4.5, 0.5, 0.5))
AlphaPlot(out3, col.gp=c('aquamarine3', 'coral3'), bgpt='beige', pch=22, legend.cex=0.8)

