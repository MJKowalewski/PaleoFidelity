## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data example-------------------------------------------------------------
library(PaleoFidelity)
str(FidData) # check the structure of the example dataset

## ----fidelity summary function------------------------------------------------
FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, report=TRUE)

## ----fidelity summary function part 2-----------------------------------------
FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat,
                report=TRUE, n.filters=30)

## ----fidelity summary function part 3-----------------------------------------
FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat,
                report=TRUE, n.filters=100)

## ----fidelity estimates-------------------------------------------------------
out1 <- FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, n.filters=30, rm.zero=F, iter=999)
out1$x
out1$y

## ----fidelity estimates outputs-----------------------------------------------
out1$xc # adjusted correlation measure summary
out1$yc # adjusted similarity measure summary

## ----fidelity estimates part 2------------------------------------------------
out1$x.stats
out1$y.stats

## ----classic fidelity plot, fig.width=7.5, fig.height=4-----------------------
par(mar=c(4, 4, 0.5, 0.5))
SJPlot(out1, gpcol=c('aquamarine3', 'coral3'), legend.cex=0.8)

## ----classic fidelity plot 2, fig.width=7.5, fig.height=4---------------------
par(mar=c(4, 4, 0.5, 0.5))
SJPlot(out1, gpcol=c('aquamarine3', 'coral3'), bubble=F, unadj=T, legend.cex=0.8)

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

## ----plot alpha 2, fig.width=7.5, fig.height=4--------------------------------
out3 <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, CI=0.95, iter=1000)
par(mar=c(4, 4.5, 0.5, 0.5))
AlphaPlot(out3, col.gp=c('aquamarine3', 'coral3'), bgpt='beige', pch=22, legend.cex=0.8)

