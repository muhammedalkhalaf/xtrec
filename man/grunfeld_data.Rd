\name{grunfeld_data}
\alias{grunfeld_data}
\title{Grunfeld Investment Data}
\description{
  A balanced panel dataset for 5 large US corporations over 20 years (1935-1954).
}
\usage{
grunfeld_data()
}
\value{
  A data frame with columns: \code{firm}, \code{year}, \code{invest}, \code{mvalue}.
}
\references{
Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment}.
PhD thesis, University of Chicago.
}
\examples{
dat <- grunfeld_data()
head(dat)
}
\keyword{datasets}
