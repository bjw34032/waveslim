\name{sine.taper}
\alias{sine.taper}
\title{Computing Sinusoidal Data Tapers}
\description{
  Computes sinusoidal data tapers directly from equations.
}
\usage{sine.taper(n, k)
}
\arguments{
  \item{n}{length of data taper(s)}
  \item{k}{number of data tapers}
}
\value{
  A vector or matrix of data tapers (cols = tapers).
}
\details{
  See reference.
}
\references{
  Riedel, K. S. and A. Sidorenko (1995)
  Minimum bias multiple taper spectral estimation,
  \emph{IEEE Transactions on Signal Processing}, \bold{43}, 188-195.
}
\seealso{
  \code{\link{dpss.taper}}.
}
%\examples{}
\author{B. Whitcher}
\keyword{ts}
