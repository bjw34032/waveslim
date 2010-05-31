\name{qmf}
\alias{qmf}
\title{Quadrature Mirror Filter}
\description{
  Computes the quadrature mirror filter from a given filter.
}
\usage{qmf(g, low2high=TRUE)
}
\arguments{
  \item{g}{Filter coefficients.}
  \item{low2high}{Logical, default is \code{TRUE} which means a low-pass
    filter is input and a high-pass filter is output.  Setting
    \code{low2high=F} performs the inverse.}
}
\value{
  Quadrature mirror filter.
}
\details{
  None.
}
\references{
  Any basic signal processing text.
}
\seealso{
  \code{\link{wave.filter}}.
}
\examples{
## Haar wavelet filter
g <- wave.filter("haar")$lpf
qmf(g)
}
\author{B. Whitcher}
\keyword{ts}
