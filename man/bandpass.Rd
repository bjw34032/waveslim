\name{Band-pass variance}
\alias{bandpass.fdp}
\alias{bandpass.spp}
\alias{bandpass.spp2}
\alias{bandpass.var.spp}
\title{Bandpass Variance for Long-Memory Processes}
\description{
  Computes the band-pass variance for fractional difference (FD) or
  seasonal persistent (SP) processes using numeric integration of their
  spectral density function.
}
\usage{
bandpass.fdp(a, b, d)
bandpass.spp(a, b, d, fG)
bandpass.spp2(a, b, d1, f1, d2, f2)
bandpass.var.spp(delta, fG, J, Basis, Length)
}
\arguments{
  \item{a}{Left-hand boundary for the definite integral.}
  \item{b}{Right-hand boundary for the definite integral.}
  \item{d,delta,d1,d2}{Fractional difference parameter.}
  \item{fG,f1,f2}{Gegenbauer frequency.}
  \item{J}{Depth of the wavelet transform.}
  \item{Basis}{Logical vector representing the adaptive basis.}
  \item{Length}{Number of elements in Basis.}
}
\value{
  Band-pass variance for the FD or SP process between \eqn{a} and
  \eqn{b}.
}
\details{
  See references.
}
\references{
  McCoy, E. J., and A. T. Walden (1996)
  Wavelet analysis and synthesis of stationary long-memory processes,
  \emph{Journal for Computational and Graphical Statistics}, \bold{5},
  No. 1, 26-56.

  Whitcher, B. (2001)
  Simulating Gaussian stationary processes with unbounded spectra,
  \emph{Journal for Computational and Graphical Statistics}, \bold{10},
  No. 1, 112-134.
}
%\seealso{}
%\examples{}
\author{Brandon Whitcher}
\keyword{ts}
