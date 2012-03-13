\name{dualfilt1}
\alias{dualfilt1}
\alias{AntonB}
\title{Kingsbury's Q-filters for the Dual-Tree Complex DWT}
\description{
  Kingsbury's Q-filters for the dual-tree complex DWT.
}
\usage{
dualfilt1()
}
\arguments{
  None.
}
\details{
  These cofficients are rounded to 8 decimal places.
}
\value{
  \item{af}{List (\eqn{i=1,2}) - analysis filters for tree \eqn{i}}
  \item{sf}{List (\eqn{i=1,2}) - synthesis filters for tree \eqn{i}}
  Note: \code{af[[2]]} is the reverse of \code{af[[1]]}.
}
\references{
  Kingsbury, N.G. (2000).  A dual-tree complex wavelet transform with
  improved orthogonality and symmetry properties, \emph{Proceedings of
    the IEEE Int. Conf. on Image Proc.} (ICIP).

  WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY\cr
  \url{http://taco.poly.edu/WaveletSoftware/}
}
\seealso{
  \code{\link{dualtree}}
}
%\examples{}
\author{Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher}
\keyword{ts}
