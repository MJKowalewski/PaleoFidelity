#' \name{FidData}
#' \docType{data}
#' \alias{FidData}
#' \title{FidData: An example of a live-dead dataset}
#' \description{
#' A live-dead dataset including live and dead specimen counts from 44 sampling station off the coast of
#' North Carolina, USA. Data include counts of live and dead macrobenthic invertebrates.
#' }
#' \usage{PaleoFidelity}
#'
#' @format A list with four objects:
#' \describe{
#'   \item{live}{matrix with 44 rows (localities) and 202 columns (species) with counts of live specimens}
#'   \item{dead}{matrix with 44 rows (localities) and 202 columns (species) with counts of dead specimens}
#'   \item{habitat}{a 2-level factor describing two main habitats}
#'   \item{fossiltype}{a 2-level factor describing categorizing species in terms their fossilization potential}
#'   ...
#' }
#' @source \url{https://www.nature.com/articles/sdata201854}
#' \keyword{datasets}
