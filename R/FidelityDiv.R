#' Live-Dead differences in alpha diversity and evenness
#'
#' FidelityDiv provides estimates of differences in diversity between pairs of live and dead samples
#' using community abundance data. In the case of datastes representing more than one site, the
#' function returns also means of differences. If 'gp' factor is provided to aggregate sets of
#' sites/samples, means for groups are returned as well.
#'
#'@details FidelityDiv assess live-dead offsets in evenness/diversity between pairs of sympatric
#' live and dead samples using a bivariate approch as described in  Olszewski and Kidwell (2007).
#' The estimates of (1) alpha diversity and (2) evenness are returned as 2 separate objects.
#'
#' (1) x - Live-dead offsets in alpha diversity for individual sites. The difference is measured
#' as the difference between natural logarithms of sample-standardized species richness of
#' sympatric dead and live samples:
#'
#' DELTA S = ln(S-dead) - log(S-live)

#' A negative value indicates that the sample standardized alpha diversity of a live sample
#' exceeds alpha diversity of the sympatric dead sample (and vice versa). Confidence intervals
#' and p.values for Null H: Delta S = 0 are also reported.
#'
#' (2) y - Live-dead offsets in evenness for individual sites. The difference is measured as the
#' difference between sample-standardized estimates of Hurlbert's PIE for live and dead samples:
#'
#' DELTA PIE = PIE(DEAD) - PIE(LIVE)
#'
#' A negative value indicates that evenness of live samples exceeds evenness of the sympatric
#' dead sample (and vice versa). Confidence intervals and p.values for Null H: Delta PIE = 0
#' are also reported.
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
#' @param report Logical (default=FALSE) to print compliance report from function FidelitySummary
#'
#' @param n.filters An integer used to filter out small samples (default n.filters=0, all samples kept)
#'
#' @param t.filters An integer used to filter out rare taxa (default t.filters=1, taxa >= 1 occurrence kept)
#'
#' @param iter An integer defining number of resampling iteration (default iter=100)
#'
#' @param CI A numerical value (default = 0.5) defining confidence bars for individual sites.
#'  Note: 0.5 - plots bars represnting inter-quartile ranges, 0.95 - plots 95% confidence bars, etc.
#'  Confidence bars are estimated as percentiles of subsampled estimates of Delta S and Delta PIE.
#'
#' @param CImean A numerical value (default = 0.99) defining confidence bars for means of all sites
#'  or groups of sites (if 'gp' factor was provided). Note: 0.5 - plots bars represnting inter-quartile
#'  range, 0.95 - plots 95% confidence bars, etc. Confidence bars are estimated as percentiles of
#'  subsampled estimates of Delta S and Delta PIE.
#'
#' @param outdata Logical (default = FALSE) to determine if data files should be included in the output
#'
#' @return A list containing the following components:
#'   \item{live}{The post-processed version of 'live' data matrix used in all analyses}
#'   \item{dead}{The post-processed version of 'dead' data matrix used in all analyses}
#'   \item{gp}{The post-processed version of 'gp' factor, when provided}
#'   \item{x}{DELTA S values for each live-dead comparisons (site-level differences in sample standardized
#'            species richness)}
#'   \item{y}{DELTA PIE values for each live-dead comparisons (site-level differences in evenness
#'            estimated as Hurlbert's PIE)}
#'   \item{xmean}{Grand mean of Delta S}
#'   \item{ymean}{Grand mean of Delta PIE}
#'   \item{xmean}{Grand mean of Delta S}
#'   \item{xgp}{Group means of Delta S (when 'gp' factor provided)}
#'   \item{ygp}{Group means of Delta PIE (when 'gp' factor provided)}
#'   \item{p.values}{p.values for Null H: Delta.S.p=0 and Delta.PIE.p=0}
#'   \item{p.gps}{per-group p.values for Null H: Delta.S.p=0 and Delta.PIE.p=0}
#'
#' @examples
#'
#' FidelityDiv(FidData$live, FidData$dead, n.filters=30)
#' my.fid <- FidelityDiv(FidData$live, FidData$dead, FidData$habitat, n.filters=50, iter=1000, CI=0.95)
#' my.fid$x # site-level estimates of Delta S with 95% CIs and p values
#' my.fid$p.gps # p values for means of groups
#'
#' @export
#'
#' @importFrom vegan rrarefy
#'
#' @references Olszewski, T.D., and Kidwell, S.M., 2007, The preservational fidelity of evenness
#'             in molluscan death assemblages. Paleobiology 33: 1:23.

FidelityDiv <- function(live, dead, gp=NULL, report=FALSE, n.filters=0, t.filters=1,
                        iter=100, CI=0.5, CImean=0.99, outdata=FALSE)
{

  # 1. Data assessment and filtering
  out <- FidelitySummary(live, dead, gp, report=report, output=TRUE,
                         n.filters=n.filters, t.filters=t.filters) # check/filter data
  if (length(out) == 2) {live <- out$live;  dead <- out$dead}
  if (length(out) == 3) {live <- out$live;  dead <- out$dead; gp <- out$gp}
  # 2. Alpha Diversity/Evenness
  min.sam <- apply(cbind(rowSums(live), rowSums(dead)), 1, min)
  pie.f <- function(x) (sum(x > 0) / (sum(x > 0) - 1)) * (1 - sum((x / sum(x)) ^ 2))
  delta.alpha <- function(x, y, min) {
    a <- vegan::rrarefy(x, sample=min)
    b <- vegan::rrarefy(y, sample=min)
    cbind(log(apply(b, 1, function(z) sum(z > 0))) - log(apply(a, 1, function(z) sum(z > 0))),
    apply(b, 1, pie.f) - apply(a, 1, pie.f))
  }
  out1 <- array(NA, dim=c(nrow(live), 2, iter))
  for (i in 1:iter) out1[,,i] <- delta.alpha(live, dead, min.sam)
  p1DS <- apply(rbind(out1[,1,]), 1, function(x) sum(x > 0))
  p2DS <- apply(rbind(out1[,1,]), 1, function(x) sum(x < 0))
  p1DP <- apply(rbind(out1[,1,]), 1, function(x) sum(x > 0))
  p2DP <- apply(rbind(out1[,1,]), 1, function(x) sum(x < 0))
  p.DS <- 2 * apply(cbind(p1DS, p2DS), 1, min) / iter
  p.DP <- 2 * apply(cbind(p1DP, p2DP), 1, min) / iter
  if (sum(p.DS == 0) > 0) p.DS[p.DS == 0] <- 1 / iter
  if (sum(p.DP == 0) > 0) p.DP[p.DP == 0] <- 1 / iter
  DS <- cbind(n.std=min.sam, est=rowMeans(rbind(out1[,1,])),
              t(apply(rbind(out1[,1,]), 1, stats::quantile, prob=c((1 - CI) / 2, 1 - (1 - CI) / 2))),
              p=p.DS)
  DP <- cbind(n.std=min.sam, est=rowMeans(rbind(out1[,2,])),
              t(apply(rbind(out1[,2,]), 1, stats::quantile, prob=c((1 - CI) / 2, 1 - (1 - CI) / 2))),
              p=p.DP)

  # 3.Means for all data and by groups (if 'gp' factor provided)
  # All data
  allS <- colMeans(rbind(out1[,1,]))
  allP <- colMeans(rbind(out1[,2,]))
  meanDS <- c(mean(allS), stats::quantile(allS, prob=c((1-CImean)/2, 1 - (1-CImean)/2)))
  meanDP <- c(mean(allP), stats::quantile(allP, prob=c((1-CImean)/2, 1 - (1-CImean)/2)))
  allpS <- 2 * min(sum(allS > 0), sum(allS < 0))
  ifelse(allpS == 0, Delta.S.p <- 1 / iter, Delta.S.p <- allpS)
  allpP <- 2 * min(sum(allP > 0), sum(allP < 0))
  ifelse(allpP == 0, Delta.PIE.p <- 1 / iter, Delta.PIE.p <- allpP)
  # By groups
  if (length(gp) > 0) {
   outDS <- sapply(as.data.frame(out1[,1,]), function(x) tapply(x, gp, mean))
   outDS2 <- cbind(n.sites=tapply(gp, gp, length), est=apply(outDS, 1, mean),
            t(apply(outDS, 1, stats::quantile, prob=c((1-CImean)/2, 1 - (1-CImean)/2))))
   outDP <- sapply(as.data.frame(out1[,2,]), function(x) tapply(x, gp, mean))
   outDP2 <- cbind(n.sites=tapply(gp, gp, length), est=apply(outDP, 1, mean),
                  t(apply(outDP, 1, stats::quantile, prob=c((1-CImean)/2, 1 - (1-CImean)/2))))
   gpDS1 <- apply(outDS, 1, function(x) sum(x<0))
   gpDS2 <- apply(outDS, 1, function(x) sum(x>0))
   gpDP1 <- apply(outDP, 1, function(x) sum(x<0))
   gpDP2 <- apply(outDP, 1, function(x) sum(x>0))
   p.gp.DS <- 2 *  apply(cbind(gpDS1, gpDS2), 1, min) / iter
   p.gp.DP <- 2 *  apply(cbind(gpDP1, gpDP2), 1, min) / iter
   if (sum(p.gp.DS == 0) > 0) p.gp.DS[p.gp.DS == 0] <- 1 / iter
   if (sum(p.gp.DP == 0) > 0) p.gp.DP[p.gp.DP == 0] <- 1 / iter
  }
outDS3 <- NULL
outDP3 <- NULL
p.GP <- NULL
if (length(gp) > 0) outDS3 <- outDS2
if (length(gp) > 0) outDP3 <- outDP2
if (length(gp) > 0)  p.GP <- cbind(p.Delta.S=p.gp.DS, p.Delta.PIE=p.gp.DP)

if (outdata) out1 <- list(live=live, dead=dead, gp=gp, out=out1, x=DS, y=DP, xmean=meanDS, ymean=meanDP,
                          xgp=outDS3, ygp=outDP3, p.values=cbind(Delta.S.p, Delta.PIE.p))
if (!outdata) out1 <- list(gp=gp, x=DS, y=DP, xmean=meanDS, ymean=meanDP, xgp=outDS3, ygp=outDP3,
                           p.values=cbind(Delta.S.p, Delta.PIE.p), p.gps=p.GP)
  return(out1)
}
