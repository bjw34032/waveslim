\name{mra.2d}
\alias{mra.2d}
\title{Multiresolution Analysis of an Image}
\description{
  This function performs a level \eqn{J} additive decomposition of the
  input matrix or image using the pyramid algorithm (Mallat 1989).
}
\usage{mra.2d(x, wf = "la8", J = 4, method = "modwt", boundary = "periodic")
}
\arguments{
  \item{x}{A matrix or image containing the data be to decomposed.  This
    must be have dyadic length in both dimensions (but not necessarily
    the same) for \code{method="dwt"}.}
  \item{wf}{Name of the wavelet filter to use in the decomposition.  By
    default this is set to \code{"la8"}, the Daubechies orthonormal
    compactly supported wavelet of length \eqn{L=8} least asymmetric
    family.}
  \item{J}{Specifies the depth of the decomposition.  This must be a
    number less than or equal to
    \eqn{\log(\mbox{length}(x),2)}{log(length(x),2)}.}
  \item{method}{Either \code{"dwt"} or \code{"modwt"}.}
  \item{boundary}{Character string specifying the boundary condition.
    If \code{boundary=="periodic"} the default, then the matrix you
    decompose is assumed to be periodic on its defined interval,\cr
    if \code{boundary=="reflection"}, the matrix beyond its boundaries
    is assumed to be a symmetric reflection of itself.}
}
\value{
  Basically, a list with the following components
  \item{LH?}{Wavelet detail image in the horizontal direction.}
  \item{HL?}{Wavelet detail image in the vertical direction.}
  \item{HH?}{Wavelet detail image in the diagonal direction.}
  \item{LL\eqn{J}}{Wavelet smooth image at the coarsest resolution.}
  \item{\eqn{J}}{Depth of the wavelet transform.}
  \item{wavelet}{Name of the wavelet filter used.}
  \item{boundary}{How the boundaries were handled.}
}
\details{
  This code implements a two-dimensional multiresolution analysis
  by performing the one-dimensional pyramid algorithm (Mallat 1989) on
  the rows and columns of the input matrix.  Either the DWT or MODWT may
  be used to compute the multiresolution analysis, which is an additive
  decomposition of the original matrix (image).
}
\references{
  Mallat, S. G. (1989)
  A theory for multiresolution signal decomposition: the wavelet
  representation,
  \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
  \bold{11}, No. 7, 674-693.

  Mallat, S. G. (1998)
  \emph{A Wavelet Tour of Signal Processing},
  Academic Press.
}
\seealso{
  \code{\link{dwt.2d}}, \code{\link{modwt.2d}}
}
\author{B. Whitcher}
\examples{
## Easy check to see if it works...
## --------------------------------

x <- matrix(rnorm(32*32), 32, 32)
# MODWT
x.mra <- mra.2d(x, method="modwt")
x.mra.sum <- x.mra[[1]]
for(j in 2:length(x.mra))
  x.mra.sum <- x.mra.sum + x.mra[[j]]
sum((x - x.mra.sum)^2)

# DWT
x.mra <- mra.2d(x, method="dwt")
x.mra.sum <- x.mra[[1]]
for(j in 2:length(x.mra))
  x.mra.sum <- x.mra.sum + x.mra[[j]]
sum((x - x.mra.sum)^2)
}
\keyword{ts}
