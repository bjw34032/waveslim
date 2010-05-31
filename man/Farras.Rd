\name{Farras}
\alias{farras}
\alias{FSfarras}
\title{Farras nearly symmetric filters}
\description{
  Farras nearly symmetric filters for orthogonal 2-channel perfect
  reconstruction filter bank and Farras filters organized for the
  dual-tree complex DWT.
}
\usage{
farras()
FSfarras()
}
\arguments{
  None.
}
\value{
  \item{af}{List ($i = 1,2$) - analysis filters for tree $i$}
  \item{sf}{List ($i = 1,2$) - synthesis filters for tree $i$}
}
\details{
}
\references{
  A. F. Abdelnour and I. W. Selesnick. "Nearly symmetric orthogonal
  wavelet bases", Proc. IEEE Int. Conf. Acoust., Speech, Signal
  Processing (ICASSP), May 2001.

  WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY\cr
  \url{http://taco.poly.edu/WaveletSoftware/}
}
\seealso{
  \code{\link{afb}}, \code{\link{dualtree}}, \code{\link{dualfilt1}}.
}
%\examples{}
\author{Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher}
\keyword{ts}
