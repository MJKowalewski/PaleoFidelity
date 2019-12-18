#' Compositional measures of live-dead fidelity
#'
#' FidelityEst estimates compositional fidelity measures by comparing two matching matrices (live
#' and dead) with community abundance data. The function returns  fidelity estimates for indiviudal
#' sites and for groups of sites if 'gp' factor is provided.
#'
#' @details FidelityEst assesses compositional fidelity using
#' measures of correlation/associations/similarity.
#' (1) x - a measure of correlation/association: spearman (default), kendall, or pearson;
#' (2) y - an abundance-based indices of similarity such as bray (default) or jaccard-chao.
#'
#' Because many of those fildeity measures tend to be sensitive to unbalanced sampling,
#' FidelityEst function attempts to correct sampling bias by assessing data-specific biases
#' in correlation/similarity measures. The bias is estimated using a resampling protocol
#' under the perfect fidelity (PF) model, in which pooled (live + dead) counts are randomly
#' partitioned into replicate pairs of samples (using sample sizes of original samples), thus
#' creating sample pairs derived from a single underlying rank abundance distribution
#' of species (i.e., perfect fidelity). For an unbiased estimator, the resampled fidelity
#' measures should indicate perfect fidelity (e.g., Spearman rho = 1). The offset between
#' the expected observed PF value (1 - PF) provides a data-specific estimate of sampling bias.
#' The adjusted fidelity measure is then given by Adjusted = Observed + (1 - PF).
#' Replicate resampling produces a distribution of PF values and resulting adjusted fidelity
#' measures, from which confidence intervals and signifiance tests can be derived.
#'
#' @param live A matrix with counts of live-collected specimens (rows=sites, columns=taxa).
#'  Dimensions of 'live' and 'dead' matrices must match exactely.
#'
#' @param dead A matrix with counts of dead-collected specimens (rows=sites, columns=taxa).
#'  Dimensions of 'live' and 'dead' matrices must match exactely.
#'
#' @param gp An optional univariate factor defining groups of sites. The length of gp must
#'  equal number of rows of 'live' and 'dead' matrices.
#'
#' @param cor.measure A character string (default='spearman') defining correlation measure
#'  (stats function 'vignettecor') used to estimate live-dead correlations.
#'
#' @param sim.measure A character string (default='chao') defining similarity measure (vegan
#'   function 'vegdist') used to estiamte live-dead similiarity. Any measure acceptable by
#'   'vegdist' can be used.
#'
#' @param n.filters An integer used to filter out small samples (default n.filters=0, all samples kept)
#'
#' @param t.filters An integer used to filter out rare taxa (default t.filters=1, taxa >= 1 occurrence kept)
#'
#' @param iter An integer defining number of resampling iteration (default iter = 99)
#'
#' @param rm.zero Logical (default rm.zero = FALSE) (removes double 0's when computing correlation measure)
#'
#' @param tfsd A character string (default='wisconsin') specifying data standardization
#' or transformations (applicable only for similarity measure).The following options are
#' available: 'none' (or any unused character string) - no standardization/transformation,
#' 'total' - relative abundance, 'wisconsin' - double relativization,
#' 'r4' - 4th root transformation, 'log' - ecological log-transformation,
#' 'total4' - 4th root transformation of relative abundances.
#'
#' @return A list containing the following components:
#'   \item{x}{Live-dead correlation coefficients for each site}
#'   \item{y}{Live-dead similarity coefficients for each live-dead comparison}
#'   \item{xc}{Adjusted live-dead correlation coefficients, 95% confidence intervals,
#'    and estimated live-dead correlation coefficients for perfect fidelity model}
#'   \item{yc}{Adjusted live-dead similarity coefficients, 95% confidence intervals,
#'    and estimated live-dead similarity coefficients for perfect fidelity model}
#'    \item{x.pf.dist}{Distributions of randomized correlation values under perfect
#'     fidelity model for each of the live-dead sample comparison}
#'    \item{y.pf.dist}{Distributions of randomized correlation values under perfect
#'     fidelity model for each of the live-dead sample comparison}
#'    \item{xc.dist}{Distributions of model adjusted correlation values for
#'    each of the samples}
#'    \item{y.pf.dist}{Distributions of model adjusted similarity values for
#'    each of the samples}
#'   \item{x.stats}{Statistical summary for correlation coefficients for all data
#'   and for each group when 'gp' factor provided}
#'   \item{y.stats}{Statistical summary for similarity coefficients for all data
#'   and for each group when 'gp' factor provided}
#'   \item{live}{The post-processed version of 'live' data matrix used in analyses}
#'   \item{dead}{The post-processed version of 'dead' data matrix used in analyses}
#'   \item{gp}{The post-processed version of 'gp', when 'gp' factor provided}
#'   \item{values}{A list with values of parameters used in the analysis}
#'
#' @examples
#'
#' data(FidData)
#' out1 <- FidelityEst(live = FidData$live[6:9,], dead = FidData$dead[6:9,],
#'                     gp = FidData$habitat[6:9], cor.measure='spearman',
#'                     sim.measure='bray', n.filters=20, iter=99, rm.zero=FALSE, tfsd='total4')
#' SJPlot(out1, gpcol=c('forestgreen', 'coral3'))
#'
#' @export
#' @importFrom vegan vegdist

FidelityEst <- function(live, dead, gp=NULL, cor.measure='spearman', sim.measure='bray',
                        n.filters=0, t.filters=1, iter=49,
                        rm.zero=FALSE, tfsd='wisconsin')
 {

# 1.1. Data assessment and filtering
  out <- FidelitySummary(live, dead, gp, report=FALSE, output=TRUE,
                         n.filters=n.filters, t.filters=t.filters) # check/filter data
  if (length(out)==2) {live <- out$live;  dead <- out$dead}
  if (length(out)==3) {live <- out$live;  dead <- out$dead; gp <- out$gp}

# 1.2. Correlation function with an option to remove double zeros
  my.cor.F <- function(x, y) {
    if(rm.zero) {
      if(sum(colSums(rbind(x, y)) == 0) > 0) {
        good.taxa <- which(colSums(rbind(x, y)) > 0)
        x <- x[good.taxa]
        y <- y[good.taxa]
        }
       }
    stats::cor(x, y, method=cor.measure)
    }

# 2. Observed fidelity values
  x1 <- as.data.frame(t(live))
  x2 <- as.data.frame(t(dead))
  cor.e <- mapply(function(x, y) my.cor.F(x, y), x1, x2)
  x3 <- x1; x4 <- x2
if(tfsd=='total') {x3 <- vegan::decostand(x1, 'total'); x4 <- vegan::decostand(x2, 'total')}
if(tfsd=='total4') {x3 <- vegan::decostand(x1, 'total')^(1/4); x4 <- vegan::decostand(x2, 'total')^(1/4)}
if(tfsd=='wisconsin') {x3 <- vegan::wisconsin(x1); x4 <- vegan::wisconsin(x2)}
if(tfsd=='log') {x3 <- vegan::decostand(x1, 'log'); x4 <- vegan::decostand(x2, 'log')}
if(tfsd=='r4') {x3 <- x1^0.25; x4 <- x2^0.25}
  sim.e <- mapply(function(x, y) 1 - vegan::vegdist(rbind(x,y), method=sim.measure), x3, x4)
  fid.sam <- data.frame(cor.e, sim.e)
  colnames(fid.sam) <- c(cor.measure, sim.measure)
  mean.measures <- c(mean(cor.e), mean(sim.e))
  names(mean.measures) <- c(cor.measure, sim.measure)
  if(length(gp) == 0)  mean.gp <- NA
  if(length(gp) > 0)  {
    mean.gp <- cbind(tapply(cor.e, gp, mean), tapply(sim.e, gp, mean))
    colnames(mean.gp) <- c(cor.measure, sim.measure)
    }

# 3. "Perfect Fidelity" MODEL
      FidPerfModel <- function(live, dead) {
        pooled <- live + dead
        rlive <- vegan::rrarefy(x=pooled, rowSums(live))
        rdead <- vegan::rrarefy(x=pooled, rowSums(dead))
        rx1 <- as.data.frame(t(rlive))
        rx2 <- as.data.frame(t(rdead))
        r.cor.e = mapply(function(x,y) my.cor.F(x, y), rx1, rx2)
        rx3 <- rx1; rx4 <- rx2
        if(tfsd=='total') {rx3 <- vegan::decostand(rx1, 'total'); rx4 <- vegan::decostand(rx2, 'total')}
        if(tfsd=='total4') {rx3 <- vegan::decostand(rx1, 'total')^(1/4); rx4 <- vegan::decostand(rx2, 'total')^(1/4)}
        if(tfsd=='wisconsin') {rx3 <- vegan::wisconsin(rx1); rx4 <- vegan::wisconsin(rx2)}
        if(tfsd=='log') {rx3 <- vegan::decostand(rx1, 'log'); rx4 <- vegan::decostand(rx2, 'log')}
        if(tfsd=='r4') {rx3 <- rx1^0.25; rx4 <- rx2^0.25}
        r.sim.e = mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), method = sim.measure), rx3, rx4)
        return(cbind(mean(r.cor.e), mean(r.sim.e)))
      }
     pf.output <- array(0, dim=c(iter, nrow(live), 2),
                        dimnames=list(1:iter, rownames(live),
                                      c(cor.measure, sim.measure)))
     for (i in 1:iter)  {
       for (j in 1:nrow(live)) {
       pf.output[i,j,] <- FidPerfModel(rbind(live[j,]), rbind(dead[j,]))
       }
     }
      perfidest <- apply(pf.output, c(2,3), stats::median)
      cor.obs.rep <- matrix(fid.sam[,1], iter, length(fid.sam[,1]), byrow=T)
      sim.obs.rep <- matrix(fid.sam[,2], iter, length(fid.sam[,2]), byrow=T)
      cor.m.adj <- (1 - pf.output[,,1]) +  cor.obs.rep
      sim.m.adj <- (1 - pf.output[,,2]) + sim.obs.rep
      cor.m.adj[cor.m.adj>1] <- 1
      sim.m.adj[sim.m.adj>1] <- 1
      cor.adj.sam <- apply(cor.m.adj, 2, mean)
      sim.adj.sam <- apply(sim.m.adj, 2, mean)
      cor.adj.sam.CI <- apply(cor.m.adj, 2, stats::quantile, prob=c(0.025, 0.975))
      sim.adj.sam.CI <- apply(sim.m.adj, 2, stats::quantile, prob=c(0.025, 0.975))
      corrected <- cbind(cor.adj.sam, sim.adj.sam)
      corrected.mean <- c(mean(cor.m.adj), mean(sim.m.adj))
      names(corrected.mean) <- c(cor.measure, sim.measure)
      x.stat <- c(mean(corrected[,1]), stats::quantile(corrected[,1], prob=c(0.025, 0.5, 0.975)))
      y.stat <- c(mean(corrected[,2]), stats::quantile(corrected[,2], prob=c(0.025, 0.5, 0.975)))
      if (length(gp) > 0) {
        statsgp.x1 <- tapply(corrected[,1], gp, mean)
        statsgp.x2 <- tapply(corrected[,1], gp, stats::quantile, prob=c(0.025, 0.5, 0.975))
        statsgp.x3 <- cbind(statsgp.x1, matrix(unlist(statsgp.x2), length(levels(gp)), 3, byrow=T))
        statsgp.y1 <- tapply(corrected[,2], gp, mean)
        statsgp.y2 <- tapply(corrected[,2], gp, stats::quantile, prob=c(0.025, 0.5, 0.975))
        statsgp.y3 <- cbind(statsgp.y1, matrix(unlist(statsgp.y2), length(levels(gp)), 3, byrow=T))
        x.stat <- rbind('all.samples'=x.stat, statsgp.x3)
        y.stat <- rbind('all.samples'=y.stat, statsgp.y3)
        colnames(x.stat) <- c('mean', '2.5%', 'median', '97.5%')
        colnames(y.stat) <- c('mean', '2.5%', 'median', '97.5%')
      }
      x.stat <- as.data.frame(x.stat)
      y.stat <- as.data.frame(y.stat)
      xc.sum <- data.frame(corrected[,1], t(cor.adj.sam.CI), perfidest[,1])
      names(xc.sum)[c(1,4)] <- c(paste(cor.measure,'.ADJ', sep=''),
                                 paste(cor.measure,'.PF', sep=''))
      yc.sum <- data.frame(corrected[,2], t(sim.adj.sam.CI), perfidest[,2])
      names(yc.sum)[c(1,4)] <- c(paste(sim.measure,'.ADJ', sep=''),
                                 paste(sim.measure,'.PF', sep=''))
### 4. OUTPUT
  out1 <- list(x=fid.sam[,1], y=fid.sam[,2], xc=xc.sum, yc=yc.sum,
               x.pf.dist = pf.output[,,1], y.pf.dist = pf.output[,,2],
               xc.dist = cor.m.adj, yc.dist = sim.m.adj, x.stats=x.stat,
               y.stats=y.stat, live = live, dead = dead, gp = gp,
               values=list(measures=c(cor.measure, sim.measure), data.transf=tfsd,
                          remove.double.zeros=rm.zero, iterations=iter))
  return(out1)
}
