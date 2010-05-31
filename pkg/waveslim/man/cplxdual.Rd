\name{Dualtree Complex}
\alias{cplxdual2D}
\alias{icplxdual2D}
\title{Dual-tree Complex 2D Discrete Wavelet Transform}
\description{
  Dual-tree complex 2D discrete wavelet transform (DWT).
}
\usage{
cplxdual2D(x, J, Faf, af)
icplxdual2D(w, J, Fsf, sf)
}
\arguments{
  \item{x}{2D array.}
  \item{w}{wavelet coefficients.}
  \item{J}{number of stages.}
  \item{Faf}{first stage analysis filters for tree \eqn{i}.}
  \item{af}{analysis filters for the remaining stages on tree \eqn{i}.}
  \item{Fsf}{last stage synthesis filters for tree \eqn{i}.}
  \item{sf}{synthesis filters for the preceeding stages.}
}
%\details{}
\value{
  For the analysis of \code{x}, the output is
  \item{w}{wavelet coefficients indexed by
    \code{[[j]][[i]][[d1]][[d2]]}, where \eqn{j=1,\ldots,J} (scale),
    \eqn{i=1} (real part) or \eqn{i=2} (imag part), \eqn{d1=1,2} and
    \eqn{d2=1,2,3} (orientations).}
  For the synthesis of \code{w}, the output is 
  \item{y}{output signal.}
}
\references{
  WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY\cr
  \url{http://taco.poly.edu/WaveletSoftware/}
}
\seealso{
  \code{\link{FSfarras}}, \code{\link{farras}}, \code{\link{afb2D}}, 
  \code{\link{sfb2D}}.
}
\examples{
\dontrun{
## EXAMPLE: cplxdual2D
x = matrix(rnorm(32*32), 32, 32)
J = 5
Faf = FSfarras()$af
Fsf = FSfarras()$sf
af = dualfilt1()$af
sf = dualfilt1()$sf
w = cplxdual2D(x, J, Faf, af)
y = icplxdual2D(w, J, Fsf, sf)
err = x - y
max(abs(err))
}
}
\author{Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher}
\keyword{ts}
