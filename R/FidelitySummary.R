#' FidelitySummary (to assess data compliance issues)
#'
#' FidelitySummary function checks if the input objects are appropriately formatted and contain
#' adequate data needed to carry out a mathematically meaningful fidelity analysis.
#'
#' @details This function is implemented in other PaleoFidelity functions. However, users are encouraged to use
#' it prior to any other functioons in order to:
#' 1). Check for errors and warnings;
#' 2). Generate a basic summary of user-provided objects (make sure to set parameter report=TRUE)
#' 3). Explore how data filtering (parameter "n.filters") affect the dimensionality of the data.
#'
#' NOTE: FidelitySummary function will not alter data objects passed through it, but only provide an initial
#' compliance evaluation. Also the function allows the user to detect small samples and assess if removing
#' those samples is advisable in terms of data dimensionality. If removal of small samples is desirable
#' for subsequent analyses, please specify the numerical value of "n.filters" in other PaleoFidelity functions.
#'
#' @param live A matrix with counts of live-collected specimens (rows=samples, columns=species/taxa)
#'
#' @param dead A matrix with counts of dead specimens (rows=samples, columns=species/taxa)
#'    (dimensions must match)
#'
#' @param gp An optional factor, with two or more levels, defining sample groups.
#'    The length of gp must equal number of rows in live and dead and at least two levels must
#'    include 2 or more observations to allow for by-group analyses.
#'
#' @param report Logical (default=FALSE), set report=TRUE to print additional warnings and data summary
#'
#' @param n.filters Numerical value used to remove small samples with n < n.filters (default value is 0)
#'
#' @return When report=TRUE, the function returns a matrix with data summary and may also return warnings.
#'     When report=FALSE, no objects are returned (only warnings may be produced)
#'
#' @examples
#' data(FidData)
#' FidelitySummary(live=FidData$live, dead=FidData$dead, report=TRUE)
#' FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, report=TRUE, n.filters=30)
#'
#' @export
#' @importFrom stats sd
#' @importFrom vegan vegdist

FidelitySummary <- function(live, dead, gp = NULL, report=FALSE, n.filters=0)

 {
  # PART 1: Initial complience checks
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
    warning('combined live + dead data contain empty columns')
  if (min(rowSums(live)) == 0)
    stop('live dataset contains empty rows')
  if (min(rowSums(dead)) == 0)
    stop('dead dataset contains empty rows')
  if (!identical(colnames(live), colnames(dead)))
    warning('column labels do not match between "live" and "dead" datasets')
  if (!identical(rownames(live), rownames(dead)))
    warning('row labels do not match between "live" and "dead" datasets')

  # Part II: Check factors
  if (length(gp) > 0)
    {
    if (length(gp) != nrow(live))
      stop('the length "gp" factor must equal the number of rows in live and dead')
    if (!is.factor(gp))
      stop('"gp" object must be a factor')
    if (sum(table(gp) == 0) > 0)
      {
      warning('empty levels detected and will be dropped')
      gp <- droplevels(gp)
      }
    if (sum(table(gp) > 1) < 2)
      warning('gp factor should include n > 1 observations for at least two levels')
    }
  if (length(gp) == 0)
    message('NOTE: gp factor has not been provided (by-group analyses and tests not possible)')

  # PART III: Apply n.filters (or not)
  if (n.filters == 0)
    message('NOTE: n.filters=0: no samples (rows) were removed')
  if (n.filters > 0)
  {
    removed <- length(unique(c(which(rowSums(live) < n.filters),
                               which(rowSums(dead) < n.filters))))
    if (removed > 0)
      {
      badsamples <-  unique(c(which(rowSums(live) < n.filters), which(rowSums(dead) < n.filters)))
      live <- live[-badsamples, ]
      dead <- dead[-badsamples, ]
      rm.samples <- length(badsamples)
      if (length(gp) > 0)   gp <- gp[-badsamples]
      }
    else
      rm.samples <- 0

    if (min(colSums(live) + colSums(dead)) == 0)
      {
      badtaxa <- which(colSums(live) + colSums(dead) == 0)
      live <- live[, -badtaxa]
      dead <- dead[, -badtaxa]
      rm.taxa <- length(badtaxa)
      }
    else
      rm.taxa <- 0

    message(paste('NOTE: n.filters=', n.filters, ', ',
      rm.samples, ' samples (rows) removed, ',
      rm.taxa, ' taxa (columns) removed',
      sep = ''))
  }

  # PART IV: Additional checks
  # check 1: check if samples with n < 30 are present in the data
  if (min(c(rowSums(live), rowSums(dead))) < 30)
    warning('small samples present (n<30): consider applying "n.filters >= 30"')
  # check 2: check again if gp factor still compliant
  if (length(gp) > 0)
  {
    if (sum(table(gp) > 1) < 2)
      warning('gp factor should include n > 1 observations for at least two levels')
  }

  # PART V: Generate report (if requested)
  if(report) # default setting "report=F" (none of the lines below executed)
  {
    ifelse(length(gp) == 0, num.groups <- 0, num.groups <- length(levels(gp)))
    ifelse(length(gp) == 0, num.use.groups <- 0, num.use.groups <- sum(table(gp)>1))
    report <- rbind('number of live samples' = nrow(live),
               'number of live taxa' = ncol(live),
               'number of dead samples' = nrow(dead),
               'number of dead taxa' = ncol(dead),
               'number of live specimens' = sum(live),
               'number of dead specimens' = sum(dead),
               'smallest sample (live)' = min(rowSums(live)),
               'smallest sample (dead)' = min(rowSums(dead)),
               'number of levels in "gp" factor' = num.groups,
                'number of observations in "gp" factor' = length(gp),
               'number of useful levels (levels with n>1)' = num.groups)
    colnames(report) <- 'outcomes'
  print(report)
  }

  if (length(gp) == 0) return(list(live=live, dead=dead))
  if (length(gp) > 0) return(list(live=live, dead=dead, gp=gp))

}
