\name{FidData}
\alias{FidData}
\docType{data}
\title{FidData: A live-dead dataset}
\description{An example of live-dead dataset based on marine macro-benthic invertebrate assemblages.}
\format{A list with four objects.
\describe{
\item{live}{matrix with 44 rows (sites) and 202 columns (species) with counts of live specimens}
\item{dead}{matrix with 44 rows (sites) and 202 columns (species) with counts of dead specimens}
\item{habitat}{a factor with two levels (habitat types)}
\item{fossiltype}{a factor with two levels (taxon type)}
 }
}
\details{A live-dead dataset including live and dead specimen counts from 44 sampling station off the coast of North Carolina, USA. Data include counts of live and dead macro-benthic invertebrates and represents a subset of a larger dataset described in Tyler and Kowalewski (2018).}
\references{Tyler, C.L., and Kowalewski, M., 2018, Regional surveys of macrobenthic shelf invertebrate communities in Onslow Bay, North Carolina, USA. Scientific Data 5, 180054.}
\source{ \url{https://www.nature.com/articles/sdata201854}}
\keyword{datasets}
\examples{str(FidData)}
\keyword{datasets}
