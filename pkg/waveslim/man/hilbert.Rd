\name{Hilbert}
\alias{dwt.hilbert}
\alias{dwt.hilbert.nondyadic}
\alias{idwt.hilbert}
\alias{modwt.hilbert}
\alias{imodwt.hilbert}
\alias{modwpt.hilbert}
\title{Discrete Hilbert Wavelet Transforms}
\description{
  The discrete Hilbert wavelet transforms (DHWTs) for seasonal and 
  time-varying time series analysis.  Transforms include the usual
  orthogonal (decimated), maximal-overlap (non-decimated) and
  maximal-overlap packet transforms.
}
\usage{
dwt.hilbert(x, wf, n.levels=4, boundary="periodic", ...)
dwt.hilbert.nondyadic(x, ...)
idwt.hilbert(y)
modwt.hilbert(x, wf, n.levels=4, boundary="periodic", ...)
imodwt.hilbert(y)
modwpt.hilbert(x, wf, n.levels=4, boundary="periodic")
}
\arguments{
  \item{x}{Real-valued time series or vector of observations.}
  \item{wf}{Hilbert wavelet pair}
  \item{n.levels}{Number of levels (depth) of the wavelet transform.}
  \item{boundary}{Boundary treatment, currently only \code{periodic} and
    \code{reflection}.}
  \item{y}{Hilbert wavelet transform object (list).}
  \item{\ldots}{Additional parametes to be passed on.}
}
%\value{}
%\details{}
\references{
  Selesnick, I. (200X). \emph{IEEE Signal Processing Magazine}

  Selesnick, I. (200X). \emph{IEEE Transactions in Signal Processing}

  Whither, B. and P.F. Craigmile (2004).  Multivariate Spectral Analysis 
  Using Hilbert Wavelet Pairs, \emph{International Journal of Wavelets, 
  Multiresolution and Information Processing}, to appear.
}
\seealso{
  \code{\link{hilbert.filter}}
}
%\examples{}
\author{B. Whitcher}
\keyword{ts}
