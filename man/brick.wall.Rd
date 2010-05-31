\name{brick.wall}
\alias{brick.wall}
\alias{dwpt.brick.wall}
\title{Replace Boundary Wavelet Coefficients with Missing Values}
\description{
  Sets the first \eqn{n} wavelet coefficients to \code{NA}.
}
\usage{brick.wall(x, wf, method="modwt")
dwpt.brick.wall(x, wf, n.levels, method="modwpt") 
}
\arguments{
  \item{x}{DWT/MODWT/DWPT/MODWPT object}
  \item{wf}{Character string; name of wavelet filter}
  \item{method}{Either \code{\link{dwt}} or \code{\link{modwt}} for
    \code{brick.wall}, or either \code{\link{dwpt}} or
    \code{\link{modwpt}} for \code{dwpt.brick.wall}}
  \item{n.levels}{depth of the wavelet transform}
}
\value{
  Same object as \code{x} only with some missing values.
}
\details{
  The fact that observed time series are finite causes boundary issues.
  One way to get around this is to simply remove any wavelet coefficient
  computed involving the boundary.  This is done here by replacing
  boundary wavelet coefficients with \code{NA}.
}
\references{
  Lindsay, R. W., D. B. Percival and D. A. Rothrock (1996).
  The discrete wavelet transform and the scale anlaysis of the
  surface properties of sea ice,
  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34},
  No.~3, 771-787.

  Percival, D. B. and A. T. Walden (2000)
  \emph{Wavelet Methods for Time Series Analysis},
  Cambridge University Press.
}
%\seealso{}
%\examples{}
\author{B. Whitcher}
\keyword{ts}
