\name{Selesnick}
\alias{cshift}
\alias{cshift2D}
\alias{pm}
\title{Miscellaneous Functions for Dual-Tree Wavelet Software}
\description{
  Miscellaneous functions for dual-tree wavelet software.
}
\usage{
cshift(x, m)
cshift2D(x, m)
pm(a, b)
}
\arguments{
  \item{x}{N-point vector}
  \item{m}{amount of shift}
  \item{a,b}{input parameters}
}
\value{
  \item{y}{vector \code{x} will be shifed by \code{m} samples to the left
    or matrix \code{x} will be shifed by \code{m} samples down.}
  \item{u}{(a + b) / sqrt(2)}
  \item{v}{(a - b) / sqrt(2)}
}
%\details{}
%\seealso{}
%\examples{}
\author{Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher}
\keyword{ts}
