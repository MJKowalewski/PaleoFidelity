#' Data filtering and compliance assessment
#'
#' FidelitySummary function evaluates if the input data are adequate and properly formatted for
#' carrying out fidelity analyses. When called by other PaleoFidelity function, FidelitySummary
#' allows the user to filter data by removing small samples and rare taxa.
#'
#' @details This function is implemented in other PaleoFidelity functions. However,
#' prior to any fidelity analysis, it is recommended  to check for errors/warnings and
#' generate a basic summary of datasets (report=TRUE), and explore how
#' filtering ("n.filters" and "t.filters") affect data dimensionality
#'
#' NOTE: FidelitySummary function provides an initial compliance evaluation and allows for
#' assessing if removing those samples is advisable. Once determined, the desired numerical
#' values of "n.filters" and "t.filters" need to be specified in other PaleoFidelity functions.
#'
#' @param live A matrix with counts of live-collected specimens (rows=sites, columns=taxa).
#'  Dimensions of 'live' and 'dead' matrices must match exactly.
#'
#' @param dead A matrix with counts of dead-collected specimens (rows=sites, columns=taxa).
#'  Dimensions of 'live' and 'dead' matrices must match exactly.
#'
#' @param gp An optional univariate factor defining groups of sites. The length of gp must
#'  equal number of rows of 'live' and 'dead' matrices.
#'
#' @param report Logical (default = FALSE), set report=TRUE to print notes,
#' warnings, and data summary
#'
#' @param n.filters Integer (default = 0) to remove small samples with n < n.filters occurrences
#'
#' @param t.filters Integer (default = 1) to remove rare taxa with t < t.filters occurrences.
#' Note that the default value of 1 keeps all taxa, but removes empty columns. This may not be
#' appropriate for certain applications such as simulations of null models or when measuring
#' live-dead rank correlation when assessing compositional fidelity.
#'
#' @param output Logical (default = FALSE) determines if an output with filtered datasets should
#' be produced.
#'
#' @return A list (returned only if output=TRUE) including the following components:
#'   \item{live}{The filtered live dataset where rows=sites and columns=taxa}
#'   \item{dead}{The filtered dead dataset where rows=sites and columns=taxa}
#'   \item{gp}{The grouping factor associated with sites (if provided)}
#'   \item{tax}{The grouping factor associated with taxa (if provided)}
#'
#' @examples
#' data(FidData)
#' FidelitySummary(live=FidData$live, dead=FidData$dead, report=TRUE)
#' FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, report=TRUE, n.filters=50)
#'
#' @export
#' @importFrom stats sd
#' @importFrom vegan vegdist

FidelitySummary <- function(live, dead, gp=NULL, report=FALSE, n.filters=0, t.filters=1,
                            output=FALSE) {

  # PART 1: Initial compliance checks
  if (sum(is.matrix(live), is.matrix(dead)) != 2)
    stop('"live" and/or "dead" object is not a matrix')
  if (sum(c(is.na(live), is.na(dead))) > 0)
    stop('missing values not allowed')
  if (!is.numeric(live) | !is.numeric(dead))
    stop('"live" and "dead" datasets must be numeric')
  if (!identical(dim(live), dim(dead)))
    stop('live and dead datasets must have the same dimensions')
  if (ncol(live) < 3)
    stop('at least 3 columns (species/variables) required to compute fidelity measures')
  if (min(colSums(live) + colSums(dead)) == 0)
  {
    if (t.filters > 0 & report == TRUE) message('NOTE: combined live+dead data contain empty columns. They will be removed')
    if (t.filters == 0 & report == TRUE) message('NOTE: combined live+dead data contain empty columns. They will NOT be removed. Use t.filters > 0')
  }
  if (min(rowSums(live)) == 0 & report == TRUE)
    warning('live dataset contains empty rows')
  if (min(rowSums(dead)) == 0 & report == TRUE)
    warning('dead dataset contains empty rows')

  if (!identical(colnames(live), colnames(dead)))
    warning('column labels do not match between "live" and "dead" datasets')
  if (!identical(rownames(live), rownames(dead)))
    warning('row labels do not match between "live" and "dead" datasets')

  # Part II: Check factors
  if (length(gp) > 0)
  {
    if (length(gp) != nrow(live))
      stop('the length of "gp" factor must equal the number of rows in live and dead')
    if (!is.factor(gp))
      stop('"gp" object must be a factor')
    if (sum(table(gp) == 0) > 0 & report == TRUE)
    {
      warning('empty levels detected and will be dropped')
      gp <- droplevels(gp)
    }
    if (sum(table(gp) > 1) < 2  & report == TRUE)
      warning('gp factor should include n > 1 observations for at least two levels')
  }
  if (length(gp) == 0  & report == TRUE)
    message('NOTE: gp factor has not been provided (by-group analyses and tests not possible)')


  # PART III: Apply n.filters and t.filters (or not)
  if (n.filters == 0 & report == TRUE)
    message('NOTE: n.filters=0: no samples (rows) were removed')

  if (t.filters == 0  & report == TRUE)
    message('NOTE: t.filters=0: no species(columns) were removed')

  if (n.filters > 0 | t.filters > 0)
  {
    removed <- length(unique(c(which(rowSums(live) < n.filters),
                               which(rowSums(dead) < n.filters))))
    if (removed == nrow(live))
      stop(paste('all samples were smaller than applied n.filters (',
                 n.filters,'): decrease n.filters', sep=''))
    if (removed > 0)
    {
      badsamples <-  unique(c(which(rowSums(live) < n.filters), which(rowSums(dead) < n.filters)))
      live <- rbind(live[-badsamples, , drop = F])
      dead <- rbind(dead[-badsamples, , drop = F])
      rm.samples <- length(badsamples)
      if (length(gp) > 0)   gp <- gp[-badsamples]
    }
    else
      rm.samples <- 0

    if (min(colSums(live>1) + colSums(dead>1)) == 0)
    {
      badtaxa <- which(colSums(live>0) + colSums(dead>0) < t.filters)
      if (length(badtaxa) > 0) {
        live <- rbind(live[, -badtaxa, drop = F])
        dead <- rbind(dead[, -badtaxa, drop = F])
        rm.taxa <- length(badtaxa)
      }
      else
        rm.taxa <- 0
    }
    else
      rm.taxa <- 0

    if (report == TRUE) { message(paste('NOTE: n.filters=', n.filters, ', ',
                                        rm.samples, ' samples (rows) removed, ',
                                        rm.taxa, ' taxa (columns) removed',
                                        sep = '')) }
  }

  # PART IV: Additional checks
  # check 1: check if samples with n < 30 are present in the data
  if (min(c(rowSums(live), rowSums(dead))) < 30 & report == T)
    warning('small samples present (n<30): consider applying "n.filters >= 30"')
  # check 2: check again if gp factor still compliant
  if (length(gp) > 0)
  {
    if (sum(table(gp) > 1) < 2 & report == TRUE)
      warning('gp factor should include n > 1 observations for at least two levels')
  }


  # PART V: Generate report (if requested)
  if(report)
  {
    ifelse(length(gp) == 0, num.groups <- 0, num.groups <- length(levels(gp)))
    ifelse(length(gp) == 0, num.use.groups <- 0, num.use.groups <- sum(table(gp)>1))
    report <- rbind('number of live samples' = nrow(live),
                    'number of dead samples' = nrow(dead),
                    'number of live taxa' = sum(colSums(live) > 0),
                    'number of dead taxa' = sum(colSums(dead) > 0),
                    'total number of taxa' = ncol(live),
                    'number of live specimens' = sum(live),
                    'number of dead specimens' = sum(dead),
                    'smallest sample (live)' = min(rowSums(live)),
                    'smallest sample (dead)' = min(rowSums(dead)),
                    'number of levels in "gp" factor' = num.groups,
                    'number of observations in "gp" factor' = length(gp),
                    'number of levels with n > 1' = num.use.groups)
    colnames(report) <- 'outcomes'
    print(report)
  }

  if(output) {
    if (length(gp) == 0) return(list(live=live, dead=dead))
    if (length(gp) > 0) return(list(live=live, dead=dead, gp=gp))
  }
}
