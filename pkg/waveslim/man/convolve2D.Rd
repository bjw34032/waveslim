\name{convolve2D}
\alias{convolve2D}
\title{Fast Column-wise Convolution of a Matrix}
\description{
  Use the Fast Fourier Transform to perform convolutions between a
  sequence and each column of a matrix.
}
\usage{
convolve2D(x, y, conj = TRUE, type = c("circular", "open"))
}
\arguments{
  \item{x}{\eqn{M{\times}N} matrix.}
  \item{y}{numeric sequence of length \eqn{N}.}
  \item{conj}{logical; if \code{TRUE}, take the complex \emph{conjugate}
     before back-transforming (default, and used for usual convolution).}
  \item{type}{character; one of \code{circular}, \code{open}
     (beginning of word is ok).  For \code{circular}, the two sequences
     are treated as \emph{circular}, i.e., periodic.

     For \code{open} and \code{filter}, the sequences are padded with
     zeros (from left and right) first; \code{filter} returns the middle 
     sub-vector of \code{open}, namely, the result of running a
     weighted mean of \code{x} with weights \code{y}.}
}
%\value{}
\details{
  This is a corrupted version of \code{convolve} made by replacing
  \code{fft} with \code{mvfft} in a few places.  It would be nice to
  submit this to the R Developers for inclusion.
}
%\references{}
\seealso{
  \code{\link{convolve}}
}
%\examples{}
\author{Brandon Whitcher}
\keyword{ts}
