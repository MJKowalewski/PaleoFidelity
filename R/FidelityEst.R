#' FidelityEst computes compositional measures of live-dead fidelity
#'
#' FidelityEst computes common live-dead fidelity measures by comparing two matching
#' matrices (live and dead) with species/taxon abundance data.
#' In the case of datastes representing more than one sample,
#' the function returns also means. If 'gp' factor is provided to define sample groups
#' means for groups are returned as well.
#'
#'@details FidelityEst assess compositional fidelity using two measures: (1) x - one of
#' the three common measures of correlation/association ('cor' in 'stats'): spearman
#' (default), kendall, or pearson; and (2) y - an abundance based measure of similarity
#' such as: chao (default) or bray ('vegdist' in 'vegan').
#'
#' @param live A species/taxon abundance matrix with counts of live-collected specimens
#'    (rows=samples, columns=species/taxa)
#'
#' @param dead A matrix with counts of dead specimens (rows=samples, columns=species/taxa)
#'    (dimensions of 'live' and 'dead' must match)
#'
#' @param gp An optional factor, with two or more levels, defining sample groups
#'    (the length of gp must equal number of rows in live and dead)
#'
#' @param cor.measure defines correlation coefficient (stats function 'vignettecor') used to measure pairwise
#'    live-dead correlations (default='spearman').
#'
#' @param sim.measure A measure of similarity to measure live-dead agreement (default='Chao')
#'
#' @param report Print compliance report from function FidelitySummary
#'
#' @param nfilters To filter out small samples. Default is nfilters=0 (all samples kept)
#'
#' @param iter Set number of iteration. The default is iter=0 (resampling simulations not carried out)
#'
#' @return A list containing the following components:
#'   \item{x}{The values of correlation coefficient for each live-dead comparison}
#'   \item{y}{The values of similarity coefficients for each live-dead comparison}
#'   \item{xlab}{The name of the correlation coefficient used}
#'   \item{ylab}{The name of the similarity coefficient used}
#'   \item{mean.x}{The mean correlation coefficient (same as 'x' if live-dead data include only one sample )}
#'   \item{mean.y}{The mean similarity coefficient (same as 'x' if live-dead data include only one sample )}
#'   \item{group.means.x}{Mean correlation coefficients by group (returned if 'gp' factor provided)}
#'   \item{group.means.y}{Mean similarity coefficients by group (returned if 'gp' factor provided)}
#'
#' @examples
#'
#' data(FidData)
#' FidelityEst(live=FidData$live, dead=FidData$dead, gp=FidData$habitat)
#'
#' @export
#' @importFrom stats cor
#' @importFrom vegan vegdist

FidelityEst <- function(live, dead, gp=NULL, cor.measure='spearman', sim.measure='chao',
                        n.filters=0, report=FALSE, iter=0, dbzero=TRUE)
 {

# 1.1. Data assessment and filtering
  out <- FidelitySummary(live, dead, gp, report=report, n.filters=n.filters) # check/filter data
  if (length(out)==2) {live <- out$live;  dead <- out$dead}
  if (length(out)==3) {live <- out$live;  dead <- out$dead; gp <- out$gp}

# 1.2. Correlation function that removes double zeros
  my.cor.F <- function(x, y, silent=TRUE)
    {
    if(dbzero)
      {
      if(sum(colSums(rbind(x, y)) == 0) > 0)
        {
        good.taxa <- which(colSums(rbind(x, y)) > 0)
        x <- x[good.taxa]
        y <- y[good.taxa]
        }
      }
    stats::cor(x,y, method=cor.measure)
    }

# 2. Observed fidelity values
  x1 <- as.data.frame(t(live))
  x2 <- as.data.frame(t(dead))
  cor.e <- mapply(function(x,y) my.cor.F(x, y), x1, x2)
  sim.e <- mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), dist=sim.measure), x1, x2)
  fid.sam <- cbind(cor.e, sim.e)
  colnames(fid.sam) <- c(cor.measure, sim.measure)
  mean.measures <- c(mean(cor.e), mean(sim.e))
  names(mean.measures) <- c(cor.measure, sim.measure)
  if(length(gp) == 0)  mean.gp <- NA
  if(length(gp) > 0)
  {
    mean.gp <- cbind(tapply(cor.e, gp, mean), tapply(sim.e, gp, mean))
    colnames(mean.gp) <- c(cor.measure, sim.measure)
  }

# 3. "Perfect Fidelity" NULL MODEL
  if (iter>0)
   {
    FidPerfModel <- function(live, dead)
     {
     pooled <- live + dead
     rlive <- vegan::rrarefy(x=pooled, rowSums(live))
     rdead <- vegan::rrarefy(x=pooled, rowSums(dead))
     rx1 <- as.data.frame(t(rlive))
     rx2 <- as.data.frame(t(rdead))
     return(cbind(r.cor.e = mapply(function(x,y) my.cor.F(x, y), rx1, rx2),
              r.sim.e = mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), dist = sim.measure), rx1, rx2)))
     }
   my.output <- array(dim=c(nrow(live), 2, iter))
   rownames(my.output) <- rownames(live)
   colnames(my.output) <- c(cor.measure, sim.measure)
   for (i in 1:iter)  my.output[,,i] <- FidPerfModel(live, dead)
   }
   ifelse(iter==0, PFMod <- NA, PFMod <- my.output)

# 4. "Zero Fidelity" NULL MODEL
   if (iter>0)
   {
     FidZeroModel <- function(live, dead)
     {
       rzlive <- t(apply(live, 1, sample))
       rzdead <- t(apply(dead, 1, sample))
       rzx1 <- as.data.frame(t(rzlive))
       rzx2 <- as.data.frame(t(rzdead))
       return(cbind(r.cor.ez = mapply(function(x,y) my.cor.F(x, y), rzx1, rzx2),
                    r.sim.ez = mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), dist = sim.measure), rzx1, rzx2)))
     }
     my.output2 <- array(dim=c(nrow(live), 2, iter))
     rownames(my.output2) <- rownames(live)
     colnames(my.output2) <- c(cor.measure, sim.measure)
     for (j in 1:iter)  my.output2[,,j] <- FidZeroModel(live,dead)
   }
   ifelse(iter==0, ZFMod <- NA, ZFMod <- my.output2)


# 5. Compute fidelity statistics
   if (iter>0)
   {
     perf.tot.means <- t(apply(PFMod, c(2,3), mean))
     perfidest <- colMeans(perf.tot.means)
     if (length(gp) > 0)
       perfidgp <- t(apply(PFMod, 3, function(x) apply(x, 2, function(y) tapply(y, gp, mean))))
     zero.tot.means <- t(apply(ZFMod, c(2,3), mean))
     zerofidest <- colMeans(zero.tot.means)
     if (length(gp) > 0)
       zerofidgp <- t(apply(ZFMod, 3, function(x) apply(x, 2, function(y) tapply(y, gp, mean))))
     if (length(gp) == 0)
       pfstatout <- list(rsam.means = perf.tot.means, perfidest = perfidest,
                       zsam.mean = zero.tot.means, zerofidest = zerofidest)
     if (length(gp) > 0)
       pfstatout <- list(rsam.means = perf.tot.means, perfidest = perfidest, perfidgp = perfidgp,
                         zsam.mean = zero.tot.means, zerofidest = zerofidest, zerofidgp = zerofidgp)
   }
   ifelse(iter==0, pfstats <- NA, pfstats <- pfstatout)


# 6. OUTPUT
   out1 <- list(sample.fid=fid.sam, observed.means=mean.measures, group.means=mean.gp,
                pfm=PFMod, zfm=ZFMod, pfstats=pfstats, live=live, dead=dead, gp=gp)
  return(out1)
}

MM1 <- ceiling(rlnorm(50, 0.1, 3))
MM2 <- MM1/sum(MM1)
out1 <- sample(1:50, replace=T, size=100, prob=MM2)
live1 <- rbind(table(factor(out1, levels=1:50)))
out2 <- sample(1:50, replace=T, size=100, prob=MM2)
dead1 <- rbind(table(factor(out2, levels=1:50)))
a <- FidelityEst(live1, dead1, report = F, iter = 100, dbzero = T, n.filters = 0)
str(a)
outt <- NULL
for (i in 1:100) {
 out1 <- sample(1:50, replace=T, size=100, prob=MM2)
 live1 <- rbind(table(factor(out1, levels=1:50)))
 out2 <- sample(1:50, replace=T, size=100, prob=MM2)
 dead1 <- rbind(table(factor(out2, levels=1:50)))
 a <- FidelityEst(live1, dead1, report=F, iter=100, dbzero=T)
 p11 <- sum(a$pfstats$rsam.means[,1] > a$sample.fid[1])
 p12 <- sum(a$pfstats$rsam.means[,1] < a$sample.fid[1])
 p21 <- sum(a$pfstats$rsam.means[,2] > a$sample.fid[2])
 p22 <- sum(a$pfstats$rsam.means[,2] < a$sample.fid[2])
 outt <- rbind(outt, cbind(2*min(p11, p12)/100, 2*min(p21, p22)/100))
# outt <- rbind(outt, cbind(2*p12/100, 2*p22/100))
}

hist(outt[,1])
hist(outt[,2])
sum(outt[,1]<=0.05)/100
sum(outt[,2]<=0.05)/100
outt
apply(outt, 2, max)
p12
a$pfstats$rsam.means[,1]
a$sample.fid[1]
plot(a$pfstats$rsam.means, xlim=c(-1,1), ylim=c(0,1), pch=16, col='blue3', cex=0.1)
 abline(v=0, h=0.5)
 points(a$pfstats$zsam.mean, col='red3', pch=16, cex=0.1)
 points(a$observed.means[1], a$observed.means[2], col='black', pch='+', cex=2)
 points(rbind(a$pfstats$perfidest), col='skyblue', pch=16, cex=0.7)
 points(rbind(a$pfstats$zerofidest), col='orange', pch=16, cex=0.7)
 points(rbind(a$pfstats$zerofidest, a$pfstats$perfidest), type='l')

 plot(a$sample.fid, xlim=c(-1,1), ylim=c(0,1))
 points(a$sample.fid, col='red')


dim(live)
l1 <- rbind(live[24,])
d1 <- rbind(dead[24,])
if(sum(colSums(rbind(l1,d1))==0)>0) {
 good.taxa <- which(colSums(rbind(l1,d1))>0)
 d1 <- rbind(d1[,good.taxa])
 l1 <- rbind(l1[,good.taxa])
}
cor(a1, a2,method='spearman')
a1 <- FidelitySummary(live, dead, gp=FidData$habitat, n.filters=30, report=TRUE)
a <- FidelityEst(live, dead, n.filters=40, gp=FidData$habitat, report=TRUE, iter=100)
a <- FidelityEst(live, dead, n.filters=10, report=TRUE, iter=10, dbzero=F)
a <- FidelityEst(l1, d1, report=TRUE, iter=100)
a$sample.fid
plot(a$sample.fid[,1], a$sample.fid[,2], xlim=c(-1,1), ylim=c(0,1))
plot(a$pfstats$rsam.means, xlim=c(-1,1), ylim=c(0,1), pch=16, col='blue3', cex=0.1)
  abline(v=0, h=0.5)
  points(a$pfstats$zsam.mean, col='red3', pch=16, cex=0.1)
  points(a$observed.means[1], a$observed.means[2], col='black', pch='+', cex=2)
  points(rbind(a$pfstats$perfidest), col='skyblue', pch=16, cex=0.7)
  points(rbind(a$pfstats$zerofidest), col='orange', pch=16, cex=0.7)
  points(rbind(a$pfstats$zerofidest, a$pfstats$perfidest), type='l')

  a$observed.means
     head(perf.tot.means)
perfidest
head(perfidgp)

t(x[,1:4])
fidmeans <- apply(a$pfm, c(1,2), mean)
fidsds <- apply(a$pfm, c(1,2), sd)
minss <- apply(cbind(rowSums(live), rowSums(dead)), 1, min)
maxss <- apply(cbind(rowSums(live), rowSums(dead)), 1, max)
mintx <- apply(cbind(rowSums(live>0), rowSums(dead>0)), 1, min)
maxtx <- apply(cbind(rowSums(live>0), rowSums(dead>0)), 1, max)
domL <- apply(live, 1, function(x) max(x)/sum(x))
domD <- apply(dead, 1, function(x) max(x)/sum(x))
domL-domD

coreval <- cbind(fidmeans, minss, maxss, maxss-minss, mintx, maxtx, maxtx-mintx, domL, domD, domL-domD)
cor(coreval)
plot(fidmeans[,2], maxss-minss, log='y')
plot(fidmeans[,2], maxtx-mintx, log='y')
plot(fidmeans[,2], domL-domD)
plot(fidmeans[,1], domL-domD)
plot(fidmeans[,1], minss, log='y')


dim(means)
a$group.means
a$x

plot(a$x, a$y, xlim=c(0.5,1), ylim=c(0.25,0.35), type='n')
 minN <- apply(cbind(rowSums(live),rowSums(dead)), 1, min)
 cexsize <- 2 * (minN - min(minN))/(max(minN)-min(minN)) + 0.5
 for (i in 1:dim(a$pfm)[3])
 {
  points(a$pfm[10,1,i], a$pfm[10,2,i], pch=21, col='skyblue', bg=adjustcolor('green', 0.3), cex=cexsize[10])
 }
 points(a$x, a$y, pch=21, cex=cexsize)
 abline(h=0.5, v=0, col='gray', lty=3)

 head(a$pfm[,,1])
hist(a$pfm, col='red')
b <- apply(a$pfm, c(1,3), mean)
dim(b)
b <- apply(a$pfm, c(2,3), mean)
dim(b)

a$pfm[,,1]
res <- apply(a$pfm, 3, mean)
dim(res)
str(a)
str(a)
help('function')
