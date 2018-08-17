#' FidelitySummary (to assess data compliance issues)
#'
#' FidelitySummary function checks if the input objects are compliant.
#' This function is used invisibly in other PaleoFidelity functions. However, users can use it directly to
#' obtain a full list of warnings and a basic summary of user-provided objects (report=T)
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
#' @param report Logical (default=F), set report=T to print additional warnings and data summary
#'
#' @return Returns errors and critical warnings.
#'     No objects returned unless report=T (return additional warnings and data summary)
#'
#' @examples
#'
#' data(FidData)
#' FidelitySummary(live=FidData$live, dead=FidData$dead, gp=FidData$habitat, report=T)



FidelitySummary <- function(live, dead, gp = NULL, report=F)
  {
  if(sum(is.matrix(live),is.matrix(dead))!=2) stop('"live" and/or "dead" object is not a matrix')

  if(sum(c(is.na(live),is.na(dead)))>0) stop('missing values not allowed')

  if (!is.numeric(live) | !is.numeric(dead))
    stop('all rows/columns in "live" and "dead" datasets must be numeric')

  if (!identical(dim(live), dim(dead)))
    stop('live and dead datasets must have the same dimensions')

  if(ncol(live)<3)
    stop('a minimum of three taxa required to compute fidelity measures')

  if (length(gp) > 0)
  {
    if (length(gp) != nrow(live))
      stop('the length "gp" factor must equal the number of rows in live and dead')
    if (!is.factor(gp))
      stop('"gp" object must be a factor')
    if (sum(table(gp)==0)>0) {
      warning('empty levels detected and will be dropped')
      gp <- droplevels(gp)
      }
    if (sum(table(gp)>1)<2)
      warning('gp factor should include n > 1 observations for at least two levels')
  }

  if (length(gp) == 0)
      warning('gp factor has not been provided')

  if (!identical(colnames(live), colnames(dead)))
      warning('column labels do not match between "live" and "dead" datasets')

  if (!identical(rownames(live), rownames(dead)))
    warning('row labels do not match between "live" and "dead" datasets')

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
               'number of groups (levels)' = num.groups,
               'number of useful groups (# levels with n > 1)' = num.groups)
  colnames(report) <- 'outcomes'
  return(report)
 }
}
