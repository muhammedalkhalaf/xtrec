\name{xtrec}
\alias{xtrec}
\title{Panel Unit Root Test Based on Recursive Detrending}
\description{
  Implements the t-REC and t-RREC panel unit root tests of Westerlund (2015).
  The t-REC assumes iid errors; the robust t-RREC accounts for serial
  correlation, cross-sectional dependence, and heteroskedasticity.
}
\usage{
xtrec(data, var, panel_id, time_id, trend = 0L,
      robust = FALSE, maxlag = -1L)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{var}{Character. Name of the variable to test.}
  \item{panel_id}{Character. Name of the panel identifier variable.}
  \item{time_id}{Character. Name of the time variable.}
  \item{trend}{Integer. Polynomial trend degree: \code{0} (constant only),
    \code{1} (linear trend), or \code{2} (quadratic trend).}
  \item{robust}{Logical. If \code{TRUE}, compute the robust t-RREC test.}
  \item{maxlag}{Integer. Maximum lag for BIC selection (-1 for automatic).}
}
\value{
  An object of class \code{"xtrec"} containing the test statistic, p-value,
  and panel summary statistics.
}
\references{
Westerlund, J. (2015). The effect of recursive detrending on panel unit root
tests. \emph{Journal of Econometrics}, 185(2), 453-467.
\doi{10.1016/j.jeconom.2014.09.013}
}
\examples{
dat <- grunfeld_data()
res <- xtrec(dat, var = "invest", panel_id = "firm",
             time_id = "year", trend = 0L, robust = FALSE)
print(res)
}
\keyword{htest}
