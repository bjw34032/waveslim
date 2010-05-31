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
  \item{x}{\eqn{N}-point vector}
  \item{m}{amount of shift}
  \item{a,b}{}
}
\value{
  \item{y}{vector \eqn{x} will be shifed by \eqn{m} samples to the left
    or matrix \eqn{x} will be shifed by \eqn{m} samples down.}
}
\details{
}
\references{
  WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY\cr
  \url{http://taco.poly.edu/WaveletSoftware/}
}
%\seealso{}
%\examples{}
\author{Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher}
\keyword{ts}
