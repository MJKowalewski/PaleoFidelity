#' FidelityEst computes compositional measures of live-dead fidelity
#'
#' FidelityEst computes common live-dead fidelity measures by comparing two matching
#' matrices (live and dead) with species/taxon abundance data.
#' In the case of datastes representing more than one sample,
#' the function returns also means. If 'gp' factor is provided to define sample groups
#' means for groups are returned as well.
#'
#'@details FidelityEst assesses compositional fidelity using three measures
#' (returned as 3 separate objects):
#'
#' (1) x - one of the three common measures of correlation/association:
#'  spearman (default), kendall, or pearson;
#'
#' (2) y - an abundance based measure of similarity such as: jaccard-chao (default)
#' or bray ('vegdist' in 'vegan');
#'
#' (3) z - Fidelity Index measured as a scaled Euclidean distance from perfect fidelity score of
#'  indices x and y, where z = 0 denotes perfect fidelity (correlation and similarity measures
#'  both equal 1); z = 1 denotes no fidelity (correlation and similarity measures both equal 0);
#'  and z < 1 denotes cases where inverse correlation is observed.
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
#' @param n.filters An integer used to filter out small samples (default n.filters=0, all samples kept)
#'
#' @param t.filters An integer used to filter out rare taxa (default t.filters=1, taxa >= 1 occurrence kept)
#'
#' @param iter An integer defining number of resampling iteration (default iter=0)
#'
#' @param dbzero Logical (default dbzero = TRUE) (removes double 0's when computing correlation measure)
#'
#' @return A list containing the following components:
#'   \item{x}{Live-dead correlation coefficients for each site}
#'   \item{y}{Live-dead similarity coefficients for each live-dead comparison}
#'   \item{z}{Live-dead fidelity index for each live-dead comparison}
#'   \item{measures}{The names of fidelity measures used}
#'   \item{observed.means}{Grand means of fidelity measures (means of  $x, $y, & $z) and
#'        by-group mean when 'gp' factor provided}
#'   \item{PF.stats}{Statistical summary under "Perfect Fidelity" model}
#'   \item{PF.dist}{Resampling distribution of fidelity measures under "Perfect Fidelity" model}
#'   \item{PF.dist.gp}{By-group resampling distribution of fidelity measures under
#'        "Perfect Fidelity" model when 'gp' factor provided}
#'   \item{gp.prob}{Pairwise significance tests for level combinations}
#'   \item{live}{The post-processed version of 'live' data matrix used in analyses}
#'   \item{dead}{The post-processed version of 'dead' data matrix used in analyses}
#'   \item{gp}{The post-processed version of 'gp', when 'gp' factor provided}
#'
#' @examples
#'
#' data(FidData)
#' out1 <- FidelityEst(live = FidData$live, dead = FidData$dead, gp = FidData$habitat)
#' SJPlot(out1, gpcol=c('forestgreen', 'coral3'))
#'
#' @export
#' @importFrom vegan vegdist

FidelityEst <- function(live, dead, gp=NULL, cor.measure='spearman', sim.measure='chao',
                        report=FALSE, n.filters=0, t.filters=1, iter=0, dbzero=TRUE)
 {

# 1.1. Data assessment and filtering
  out <- FidelitySummary(live, dead, gp, report=report, n.filters=n.filters, t.filters=t.filters) # check/filter data
  if (length(out)==2) {live <- out$live;  dead <- out$dead}
  if (length(out)==3) {live <- out$live;  dead <- out$dead; gp <- out$gp}

# 1.2. Correlation function with an option to remove double zeros
  my.cor.F <- function(x, y) {
    if(dbzero) {
      if(sum(colSums(rbind(x, y)) == 0) > 0) {
        good.taxa <- which(colSums(rbind(x, y)) > 0)
        x <- x[good.taxa]
        y <- y[good.taxa]
        }
       }
    stats::cor(x, y, method=cor.measure)
    }

# 2. Observed fidelity values
  fid.M <- function(x, y) sqrt((1 - x)^2 + (1 - y)^2) / sqrt(2) # fidelity index function
  x1 <- as.data.frame(t(live))
  x2 <- as.data.frame(t(dead))
  cor.e <- mapply(function(x, y) my.cor.F(x, y), x1, x2)
  sim.e <- mapply(function(x, y) 1 - vegan::vegdist(rbind(x,y), method=sim.measure), x1, x2)
  fid.e <- fid.M(cor.e, sim.e)
  fid.sam <- cbind(cor.e, sim.e, fid.e)
  colnames(fid.sam) <- c(cor.measure, sim.measure, "fid.index")
  mean.measures <- c(mean(cor.e), mean(sim.e), mean(fid.e))
  names(mean.measures) <- c(cor.measure, sim.measure, "fid.index")
  if(length(gp) == 0)  mean.gp <- NA
  if(length(gp) > 0)  {
    mean.gp <- cbind(tapply(cor.e, gp, mean), tapply(sim.e, gp, mean), tapply(fid.e, gp, mean))
    colnames(mean.gp) <- c(cor.measure, sim.measure, "fid.index")
    }

# 3. "Perfect Fidelity" NULL MODEL
    if (iter > 0) {
      FidPerfModel <- function(live, dead) {
        pooled <- live + dead
        rlive <- vegan::rrarefy(x=pooled, rowSums(live))
        rdead <- vegan::rrarefy(x=pooled, rowSums(dead))
        rx1 <- as.data.frame(t(rlive))
        rx2 <- as.data.frame(t(rdead))
        r.cor.e = mapply(function(x,y) my.cor.F(x, y), rx1, rx2)
        r.sim.e = mapply(function(x,y) 1-vegan::vegdist(rbind(x,y), method = sim.measure), rx1, rx2)
        return(cbind(mean(r.cor.e), mean(r.sim.e), mean(fid.M(r.cor.e, r.sim.e))))
      }
     my.output <- matrix(0, iter, 3)
     for (i in 1:iter)  my.output[i,] <- FidPerfModel(live, dead)
     colnames(my.output) <- c(cor.measure, sim.measure, 'Fid.index')
     perfidest <- colMeans(my.output)
     ppfFI.1 <- sum(my.output[,3] >= mean.measures[3])
     ppfFI.2 <- sum(my.output[,3] <= mean.measures[3])
     pcorr <- (2 * min(ppfFI.1, ppfFI.2) + 1) / (iter + 1)
     pfstatout <- data.frame(group = 'all data', num.sites = rep(nrow(live), 3), PFest = perfidest,
                       observed = mean.measures, PF.debt = perfidest - mean.measures,
                       p = c(NA, NA, pcorr))
    }
  ifelse(exists('my.output'), PF.stats <- pfstatout, PF.stats <- NA)
  ifelse(exists('my.output'), PF.dist <- my.output, PF.dist <- NA)

# 4. STATISTICAL COMPARISON OF GROUPS (PAIRWISE COMPARISONS)
  if (length(gp)>0 & iter > 0) {
    outlist2 <- vector("list", length(levels(gp)))
    kk <- 0
    for (k in levels(gp)) {
    kk <- kk + 1
    my.output.gp <- matrix(0, iter, 3)
    colnames(my.output.gp) <- c(cor.measure, sim.measure, 'Fid.index')
    livegp <- live[which(gp == k),]
    deadgp <- dead[which(gp == k),]
    for (j in 1:iter)  my.output.gp[j,] <- FidPerfModel(livegp, deadgp)
    perfidest.gp <- colMeans(my.output.gp)
    ppfFI.gp.1 <- sum(my.output.gp[,3] >= mean.gp[kk, 3])
    ppfFI.gp.2 <- sum(my.output.gp[,3] <= mean.gp[kk, 3])
    pcorr.gp <- (1 + 2 * min(ppfFI.gp.1, ppfFI.gp.2)) / (iter + 1)
    pfstatout.gp <- data.frame(group = k, num.sites = rep(table(gp)[k], 3),
                               PFest = perfidest.gp, observed = mean.gp[kk,],
                               PF.debt = perfidest.gp - mean.gp[kk,], p = c(NA, NA, pcorr.gp))
    ifelse(kk == 1, outlist1 <- pfstatout.gp, outlist1 <- rbind(outlist1, pfstatout.gp))
    outlist2[[kk]] <- my.output.gp
    }
    p.gp.out <- NULL
    for (m in 1:length(levels(gp))) {
      for (n in 2:length(levels(gp))) {
         if (n > m) {
      bgp1 <- outlist2[[m]][,3] - mean.gp[m,3]
      bgp2 <- outlist2[[n]][,3] - mean.gp[n,3]
      p.gp <- (1 + (2 * min(sum(bgp1 >= bgp2), sum(bgp1 <= bgp2))))/(iter+1)
      num.comp <- (length(levels(gp))/2) * (length(levels(gp)) - 1)
      ifelse(p.gp < 0.05, sgnf.1 <- '***', sgnf.1 <- '-')
      ifelse(p.gp*num.comp < 0.05, sgnf.2 <- '***', sgnf.2 <- '-')
      p.gp.out <- rbind(p.gp.out, data.frame(gp1=levels(gp)[m], gp2=levels(gp)[n],
                        p.value=p.gp, sign=sgnf.1, Bonferroni=sgnf.2))
      }
     }
    }
   }
  ifelse(exists('my.output.gp'), PF.stats.gp <- outlist1, PF.stats.gp <- NA)
  ifelse(exists('my.output.gp'), PF.dist.gp <- outlist2, PF.dist.gp <- NA)
  ifelse(exists('my.output.gp'), PF.gp.p <- p.gp.out, PF.gp.p <- NA)

# 5. OUTPUT
  ifelse(length(gp) > 0, rep1 <- rbind(total=mean.measures, mean.gp), rep1 <- mean.measures)
  ifelse(length(gp) > 0, rep2 <- rbind(PF.stats, PF.stats.gp), rep2 <- PF.stats)
  out1 <- list(x=fid.sam[,1], y=fid.sam[,2], z=fid.sam[,3],
                measures=c(cor.measure, sim.measure, "fid.index"),
                observed.means = rep1, PF.stats = rep2, PF.dist = PF.dist,
                PF.dist.gp = PF.dist.gp, gp.prob = PF.gp.p,
                live = live, dead = dead, gp = gp)
  return(out1)
}
