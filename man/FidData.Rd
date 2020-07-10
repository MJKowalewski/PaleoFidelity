\name{FidData}
\alias{FidData}
\docType{data}
\title{
FidData: An example of a live-dead dataset
}
\description{
A live-dead dataset including live and dead specimen counts from 44 sampling station off
the coast of North Carolina, USA.
}
\usage{data("FidData")}
\format{
  A list with the following 4 objects.
  \describe{
 \item{live}{matrix with 44 rows (localities) and 202 columns (species) with counts of live specimens}
 \item{dead}{matrix with 44 rows (localities) and 202 columns (species) with counts of dead specimens}
 \item{habitat}{a 2-level factor describing two main habitats}
 \item{fossiltype}{a 2-level factor categorizing species in terms of their fossilization potential}
 }
}
\details{
A live-dead dataset including live and dead specimen counts from 44 sampling station off
the coast of North Carolina, USA. Data include counts of live and dead macrobenthic invertebrates.
For more information about this dataset please see Tyler and Kowalewski (2018).
}
\source{
\url{https://www.nature.com/articles/sdata201854}
}
\references{
Tyler, C.L., Kowalewski, M., 2018, Regional surveys of macrobenthic shelf invertebrate communities in Onslow Bay, North Carolina, U.S.A. Scientific Data 5:180054 doi: 10.1038/sdata.2018.54.
}
\examples{
data(FidData)
str(FidData)
}
\keyword{datasets}
