\name{hdr}
\alias{hdr}
\alias{hdr.boxplot}
\alias{hdr.den}
\title{Highest Density Regions}
\description{
Calculates and plots highest density regions in one dimension including the HDR boxplot.
}
\usage{
hdr(x, prob = c(50, 95, 99), den, h, nn=5000, all.modes=FALSE)
hdr.den(x, prob = c(50, 95, 99), den, h, xlab="", ylab="", ...)
hdr.boxplot(x, prob = c(99, 50), h, boxlabels = "", col = gray((9:1)/10), 
    main="", xlab="", ylab="", pch=1, ...)
}
\arguments{
  \item{x}{Numeric vector containing data. In \code{hdr} and \code{hdr.den},
    if \code{x} is missing then \code{den} must be provided, and the
    HDR is computed from the given density.
    For \code{hdr.boxplot}, \code{x} can be a list containing several vectors.}
  \item{prob}{Probability coverage required for HDRs}
  \item{den}{Density of data as list with components \code{x} and \code{y}.
    If omitted, the density is estimated from \code{x} using \code{\link[base]{density}}.}
  \item{h}{Optional bandwidth for calculation of density.}
  \item{nn}{Number of random numbers used in computing f-alpha quantiles.}
  \item{all.modes}{Return all local modes or just the global mode?}
  \item{boxlabels}{Label for each box plotted.}
  \item{col}{Colours for regions of each box.}
  \item{main}{Overall title for the plot.}
  \item{xlab}{Label for x-axis.}
  \item{ylab}{Label for y-axis.}
  \item{pch}{Plotting character.}
  \item{\dots}{Other arguments passed to plot.}
}
\details{Either \code{x} or \code{den} must be provided. Hyndman's (1996) density quantile
algorithm is used for calculation. \code{hdr.den} plots
the density with the HDRs superimposed. \code{hdr.boxplot} displays a boxplot based on HDRs.
}
\value{\code{hdr.boxplot} retuns nothing. \code{hdr} and \code{hdr.den} return a list of three components:
  \item{hdr}{The endpoints of each interval in each HDR}
  \item{mode}{The estimated mode of the density.}
  \item{falpha}{The value of the density at the boundaries of each HDR.}
}
\references{Hyndman, R.J. (1996) Computing and graphing highest density regions
\emph{American Statistician}, \bold{50}, 120-126.}

\author{Rob Hyndman}

\seealso{\code{\link{hdr.boxplot.2d}}}
\examples{
hdr(faithful$eruptions)
hdr.boxplot(faithful$eruptions)
hdr.den(faithful$eruptions)

x <- c(rnorm(100,0,1), rnorm(100,3,1))
par(mfrow=c(1,2))
boxplot(x)
hdr.boxplot(x)
par(mfrow=c(1,1))
hdr.den(x)
}
\keyword{smooth}
\keyword{distribution}
\keyword{hplot}
