% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dwpt.R
\name{css.test}
\alias{css.test}
\alias{cpgram.test}
\alias{entropy.test}
\alias{portmanteau.test}
\title{Testing the Wavelet Packet Tree for White Noise}
\usage{
cpgram.test(y, p = 0.05, taper = 0.1)

css.test(y)

entropy.test(y)

portmanteau.test(y, p = 0.05, type = "Box-Pierce")
}
\arguments{
\item{y}{wavelet packet tree (from the DWPT)}

\item{p}{significance level}

\item{taper}{weight of cosine bell taper (\code{cpgram.test} only)}

\item{type}{\code{"Box-Pierce"} and \code{other} recognized
(\code{portmanteau.test} only)}
}
\value{
Boolean vector of the same length as the number of nodes in the
wavelet packet tree.
}
\description{
A wavelet packet tree, from the discrete wavelet packet transform (DWPT), is
tested node-by-node for white noise.  This is the first step in selecting an
orthonormal basis for the DWPT.
}
\details{
Top-down recursive testing of the wavelet packet tree is
}
\examples{

data(mexm)
J <- 6
wf <- "la8"
mexm.dwpt <- dwpt(mexm[-c(1:4)], wf, J)
## Not implemented yet
## plot.dwpt(x.dwpt, J)
mexm.dwpt.bw <- dwpt.brick.wall(mexm.dwpt, wf, 6, method="dwpt")
mexm.tree <- ortho.basis(portmanteau.test(mexm.dwpt.bw, p=0.025))
## Not implemented yet
## plot.basis(mexm.tree)

}
\references{
Brockwell and Davis (1991) \emph{Time Series: Theory and
Methods}, (2nd. edition), Springer-Verlag.

Brown, Durbin and Evans (1975) Techniques for testing the constancy of
regression relationships over time, \emph{Journal of the Royal Statistical
Society B}, \bold{37}, 149-163.

Percival, D. B., and A. T. Walden (1993) \emph{Spectral Analysis for
Physical Applications: Multitaper and Conventional Univariate Techniques},
Cambridge University Press.
}
\seealso{
\code{\link{ortho.basis}}.
}
\author{
B. Whitcher
}
\keyword{ts}
