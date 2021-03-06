% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdr.R
\name{hdr}
\alias{hdr}
\title{Highest Density Regions}
\usage{
hdr(
  x = NULL,
  prob = c(50, 95, 99),
  den = NULL,
  h = hdrbw(BoxCox(x, lambda), mean(prob)),
  lambda = 1,
  nn = 5000,
  all.modes = FALSE
)
}
\arguments{
\item{x}{Numeric vector containing data. If \code{x} is missing then
\code{den} must be provided, and the HDR is computed from the given density.}

\item{prob}{Probability coverage required for HDRs}

\item{den}{Density of data as list with components \code{x} and \code{y}.
If omitted, the density is estimated from \code{x} using
\code{\link[stats]{density}}.}

\item{h}{Optional bandwidth for calculation of density.}

\item{lambda}{Box-Cox transformation parameter where \code{0 <= lambda <=
1}.}

\item{nn}{Number of random numbers used in computing f-alpha quantiles.}

\item{all.modes}{Return all local modes or just the global mode?}
}
\value{
A list of three components: \item{hdr}{The endpoints of each interval
in each HDR} \item{mode}{The estimated mode of the density.}
\item{falpha}{The value of the density at the boundaries of each HDR.}
}
\description{
Calculates highest density regions in one dimension
}
\details{
Either \code{x} or \code{den} must be provided. When \code{x} is provided,
the density is estimated using kernel density estimation. A Box-Cox
transformation is used if \code{lambda!=1}, as described in Wand, Marron and
Ruppert (1991). This allows the density estimate to be non-zero only on the
positive real line. The default kernel bandwidth \code{h} is selected using
the algorithm of Samworth and Wand (2010).

Hyndman's (1996) density quantile algorithm is used for calculation.
}
\examples{
# Old faithful eruption duration times
hdr(faithful$eruptions)
}
\references{
Hyndman, R.J. (1996) Computing and graphing highest density
regions. \emph{American Statistician}, \bold{50}, 120-126.

Samworth, R.J. and Wand, M.P. (2010). Asymptotics and optimal bandwidth
selection for highest density region estimation.  \emph{The Annals of
Statistics}, \bold{38}, 1767-1792.

Wand, M.P., Marron, J S., Ruppert, D. (1991) Transformations in density
estimation. \emph{Journal of the American Statistical Association},
\bold{86}, 343-353.
}
\seealso{
\code{\link{hdr.den}}, \code{\link{hdr.boxplot}}
}
\author{
Rob J Hyndman
}
\keyword{distribution}
\keyword{smooth}
