## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data example--------------------------------------------------------
library(PaleoFidelity)
str(FidData) # check the structure of the example dataset

## ----fidelity estimates--------------------------------------------------
out1 <- FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, n.filters=30, dbzero=F, iter=999)
out1$observed.means

## ----fidelity estimates part 2-------------------------------------------
out1$PF.stats

## ----fidelity estimates part 3-------------------------------------------
out1$gp.prob

## ----classic fidelity plot-----------------------------------------------
SJPlot(out1, gpcol=c('aquamarine3', 'coral3'))

## ----classic fidelity plot 2---------------------------------------------
SJPlot(out1, gpcol=c('aquamarine3', 'coral3'), bubble=F)

## ----alpha 1-------------------------------------------------------------
out2 <- FidelityDiv(FidData$live, FidData$dead, iter=1000)
str(out2)

## ----plot alpha 1--------------------------------------------------------
AlphaPlot(out2)

## ----alpha 2-------------------------------------------------------------
out3 <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, CI=0.95, iter=1000)
str(out3)

## ----plot alpha 2--------------------------------------------------------
AlphaPlot(out3, col.gp=c('aquamarine3', 'coral3'), bgpt='beige', pch=22)

