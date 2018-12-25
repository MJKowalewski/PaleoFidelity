#' Fidelity Data example
#'
#' \name{FidData}
#' \docType{data}
#' \alias{FidData}
#' \title{FidData: An example of a live-dead dataset}
#' \description{
#' An example of live-dead dataset based on marine macrobenthic invertebrate assemblages.
#' }
#' \details{
#' A live-dead dataset including live and dead specimen counts from 44 sampling station off
#' the coast of North Carolina, USA. Data include counts of live and dead macrobenthic invertebrates
#' and represents a part of the dataset described in Tyler and Kowalewski (2018).
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
#' @references Tyler, C.L., and Kowalewski, M., 2018, Regional surveys of macrobenthic shelf
#'  invertebrate communities in Onslow Bay, North Carolina, USA. Scientific Data 5, 180054.
#' @source \url{https://www.nature.com/articles/sdata201854}
#' \keyword{datasets}
