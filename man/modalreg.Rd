\name{modalreg}
\alias{modalreg}

\title{Nonparametric Multimodal  Regression}
\description{
  Nonparametric multi-valued regression based on the modes of conditional density estimates.
}
\usage{
modalreg(x, y, xfix=seq(min(x),max(x),l=50), a, b, deg = 0, iter = 30, P = 2,
    start = "e", shrink = TRUE, shrink.const = 10, plot.type = c("p", 1),
    labels = c("", "x", "y"), pch=".", ...)
}

\arguments{
  \item{x}{Numerical vector: the conditioning variable.}
  \item{y}{Numerical vector: the response variable.}
  \item{xfix}{Numerical vector corresponding to the input values of which the fitted values shall be calculated.}
  \item{a}{Optional bandwidth in \eqn{x}-direction}
  \item{b}{Optional bandwidth in \eqn{y}-direction }
  \item{deg}{Degree of local polynomial used in estimation (0 or 1).}
  \item{iter}{Positive integer giving the number of mean shift iterations per point and branch.}
  \item{P}{ Maximal number of branches.}
  \item{start}{Character determining how the starting points are selected.
  \code{"q"}: proportional to quantiles; \code{"e"}: equidistant; \code{"r"}: random.
   All, \code{"q"}, \code{"e"}, and \code{"r"}, give constant starting points.
   As an alternative, the choice \code{"v"} gives variable starting points equal
    to \code{"q"} for the smallest \code{x}, and equal to
   the previously fitted values for all subsequent \code{x}. }
  \item{shrink}{ Boolean. If TRUE, parts of branches are dismissed where their
  associated kernel density value falls below the threshold
   \code{1/(shrink.const*(max(x)-min(x))*(max(y)-min(y)))} }
  \item{shrink.const}{ Numerical value giving the constant used above (the higher, the less shrinking)}
  \item{plot.type}{Vector with two elements. The first one is character-valued,
   with possible values "p",  "l", and "n".  If equal to "n", no plotted output
   is given at all. If equal to "p", fitted curves are symbolized as points in the
  graphical output, otherwise as lines. (Warning: The latter setting might
  give confusing outputs especially for high \code{P}.)}.  The  second vector
  component is a numerical value either being \eqn{0} or \eqn {1}. If \eqn{1},
  the position of the starting point is depicted in the plot, otherwise omitted.
  \item{labels}{Vector of three character strings.
  The first one is the "main" title of the graphical output,
  the second one is the label of the \eqn{x} axis, and the third one the label of the
   \eqn{y} axis.}
  \item{pch}{Plotting character.}
  \item{\dots}{ Other arguments passed to \code{\link{cde.bandwidths}}.}
}
\details{
    Computes multi-modal nonparametric regression curves based on the
    maxima of conditional density estimates. The tool for the estimation
    is the conditional mean shift as outlined in Einbeck and Tutz (2006).
    Estimates of the conditional modes might fluctuate highly if \code{deg=1}.
     Hence, \code{deg=0}  is strongly recommended. For bandwidth selection, the
      hybrid rule introduced by Bashtannyk and Hyndman (2001) is employed
      if \code{deg=0}; and the Hyndman and Yao (2002) algorithm if \code{deg=1}.


}
\value{
A list with the following components:
  \item{Fitted.values}{A \code{[P x length(xfix)]}- matrix with fitted j-th branch
  in the j-th row (\eqn{1 \le j \le  P}{1 <=j <=P})
  }
  \item{Bandwidths}{A vector with bandwidths \code{a} and \code{b}}
  \item{Density}{A \code{[P x length(xfix)]}- matrix with estimated kernel densities}
  \item{Threshold}{The shrinkage threshold.}
}

\references{
 Einbeck, J., and Tutz, G. (2006) "Modelling beyond regression functions:
 an application of multimodal regression to speed-flow data". Applied Statistics,
 2006, in press.

Bashtannyk, D.M., and Hyndman, R.J. (2001) "Bandwidth selection for kernel
conditional density estimation". \emph{Computational Statistics and Data Analysis},
 \bold{36}(3), 279-298.

Hyndman, R.J. and Yao, Q. (2002) "Nonparametric estimation and
symmetry tests for conditional density functions". \emph{Journal of
Nonparametric Statistics}, \bold{14}(3), 259-278.

 }

\author{
  Jochen Einbeck (2006)
}

\seealso{ \code{\link{cde.bandwidths}}}

\examples{
  lane2.fit <- modalreg(lane2$flow, lane2$speed, xfix=(1:55)*40, a=100, b=4)
}

\keyword{regression}
\keyword{nonparametric}
