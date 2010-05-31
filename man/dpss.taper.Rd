\name{dpss.taper}
\alias{dpss.taper}
\title{Calculating Thomson's Spectral Multitapers by Inverse Iteration}
\description{
  The following function links the subroutines in "bell-p-w.o" to an 
  R function in order to compute discrete prolate spheroidal 
  sequences (dpss).
}
\usage{dpss.taper(n, k, nw = 4, nmax = 2^(ceiling(log(n, 2))))
}
\arguments{
  \item{n}{length of data taper(s)}
  \item{k}{number of data tapers; 1, 2, 3, ... (do not use 0!)}
  \item{nw}{product of length and half-bandwidth parameter (w)}
  \item{nmax}{maximum possible taper length, necessary for FORTRAN code}
}
\value{
  \item{v}{matrix of data tapers (cols = tapers)}
  \item{eigen}{eigenvalue associated with each data taper}
  \item{iter}{total number of iterations performed}
  \item{n}{same as input}
  \item{w}{half-bandwidth parameter}
  \item{ifault}{0 indicates success, see documentation for "bell-p-w"
    for information on non-zero values}
}
\details{
  Spectral estimation using a set of orthogonal tapers is becoming
  widely used and appreciated in scientific research.  It produces direct
  spectral estimates with more than 2 df at each Fourier frequency,
  resulting in spectral estimators with reduced variance.  Computation of
  the orthogonal tapers from the basic defining equation is difficult,
  however, due to the instability of the calculations -- the eigenproblem
  is very poorly conditioned.  In this article the severe numerical
  instability problems are illustrated and then a technique for stable
  calculation of the tapers -- namely, inverse iteration -- is
  described. Each iteration involves the solution of a matrix equation.
  Because the matrix has Toeplitz form, the Levinson recursions are used
  to rapidly solve the matrix equation.  FORTRAN code for this method is
  available through the Statlib archive.  An alternative stable method is
  also briefly reviewed.
}
\references{
  B. Bell, D. B. Percival, and A. T. Walden (1993)
  Calculating Thomson's spectral multitapers by inverse iteration,
  \emph{Journal of Computational and Graphical Statistics},
  \bold{2}, No. 1, 119-130.

  Percival, D. B. and A. T. Walden (1993)
  \emph{Spectral Estimation for Physical Applications: Multitaper and
    Conventional Univariate Techniques},
  Cambridge University Press.

}
\seealso{
  \code{\link{sine.taper}}.
}
%\examples{}
\author{B. Whitcher}
\keyword{ts}
