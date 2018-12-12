#' FidelityDiv estimates live-dead differences in alpha and beta diversity
#'
#' FidelityDiv provides estimates of differences in diversity between live and dead samples
#' using matrices (live and dead) with species/taxon abundance data. In the case of datastes
#' representing more than one sample, the function returns also means of differences.
#' If 'gp' factor is provided to aggregate sets of sites/samples, means for groups are returned as well.
#'
#'@details FidelityDiv assess live-dead offsets in evenness/diversity using measures
#' of alpha diversity, alpha evenness, and beta diversity (returned as 3 separate objects):
#'
#' (1) x - Live-dead offsets in alpha diversity for individual sites. The difference is measured as
#' ln(S)DEAD - ln(S)LIVE (i.e., difference between natural logarithms of sample-standardized
#' species richness of dead and live samples). A negative value indicates that alpha diversity
#' of live samples exceeds alpha diversity of dead sample (and vice versa).
#'
#' (2) y - Live-dead offsets in evenness for individual sites. The difference is measured as
#' PIE(DEAD) - PIE(LIVE) (i.e., difference between estimates of Hurlbert's PIE for live and
#' dead samples). A negative value indicates that evenness of live samples exceeds evenness
#' of dead sample (and vice versa).
#'
#' (3) z - Live-dead offsets in beta diversity for individual sites. The offset is estimated as
#'  Beta(DEAD) - Beta(LIVE). Three  measures of beta diversity can be selected:
#'  "beta variance", "Shannon's beta", "Whittaker's beta". A negative value indicates that
#'  beta diversity of live sample exceed beta diversity of dead sample (and vice versa).
#'
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
#' @param report Logical (default=FALSE) to print compliance report from function FidelitySummary
#'
#' @param n.filters An integer used to filter out small samples (default n.filters=0, all samples kept)
#'
#' @param t.filters An integer used to filter out rare taxa (default t.filters=1, taxa >= 1 occurrence kept)
#'
#' @param iter An integer defining number of resampling iteration (default iter=0)
#'
#'
#' @return A list containing the following components:
#'   \item{x}{DELTA S values for each live-dead comparisons (site-level differences in sample standardized
#'            species richness)}
#'   \item{y}{DELTA PIE values for each live-dead comparisons (site-level differences in evenness
#'            estimated as Hurlbert's PIE)}
#'   \item{z}{DELTA BETA Live-dead fidelity index for each live-dead comparison}
#'   \item{observed.means}{Grand means of diversity measures (means of  $x, $y, & $z) and
#'        by-group means when gp factor provided}
#'   \item{live}{The post-processed version of 'live' data matrix used in all analyses}
#'   \item{dead}{The post-processed version of 'dead' data matrix used in all analyses}
#'   \item{gp}{The post-processed version of 'gp' factor, when provided}
#'
#' @examples
#'
#' data(FidData)
#' out1 <- FidelityDiv(live = FidData$live, dead = FidData$dead, gp = FidData$habitat)
#'
#'
#' @export
#'
#' @importFrom vegan rrarefy

FidelityDiv <- function(live, dead, gp=NULL, report=FALSE, n.filters=0, t.filters=1, iter=100)
{

  # 1. Data assessment and filtering
  out <- FidelitySummary(live, dead, gp, report=report, n.filters=n.filters, t.filters=t.filters) # check/filter data
  if (length(out)==2) {live <- out$live;  dead <- out$dead}
  if (length(out)==3) {live <- out$live;  dead <- out$dead; gp <- out$gp}

  # 2. Alpha Diversity
  min.sam <- apply(cbind(rowSums(live), rowSums(dead)), 1, min)
  pie.f <- function(x) (sum(x>0) / (sum(x>0) - 1)) * (1 - sum((x/sum(x))^2))
  delta.alpha <- function(x, y, min) {
    a <- vegan::rrarefy(x, sample=min)
    b <- vegan::rrarefy(y, sample=min)
    cbind(log(apply(b, 1, function(z) sum(z > 0))) - log(apply(a, 1, function(z) sum(z > 0))),
    apply(b, 1, pie.f) - apply(a, 1, pie.f))
  }
  out.alpha <- array(NA, dim=c(nrow(live), 2, iter))
  for (i in 1:iter) out.alpha[,,i] <- delta.alpha(live, dead, min.sam)

  # 3. OUTPUT
  #  ifelse(length(gp) > 0, rep1 <- rbind(total=mean.measures, mean.gp), rep1 <- mean.measures)
  #  ifelse(length(gp) > 0, rep2 <- rbind(PF.stats, PF.stats.gp), rep2 <- PF.stats)
  out1 <- list(x='x', y='y', z='z', observed.means='obs.mean', live=live, dead=dead, gp=gp,
               out=out.alpha)
  return(out1)
}

