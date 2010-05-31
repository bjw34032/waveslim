\name{shift.2d}
\alias{shift.2d}
\title{Circularly Shift Matrices from a 2D MODWT}
\description{
  Compute phase shifts for wavelet sub-matrices based on the ``center of
  energy'' argument of Hess-Nielsen and Wickerhauser (1996).
}
\usage{shift.2d(z, wf, inv=FALSE)
}
\arguments{
  \item{z}{2D MODWT object}
  \item{wf}{Character string for wavelet filter.}
  \item{inv}{Boolean value on whether to perform the forward or inverse
    operation.}
}
\value{
  2D MODWT object with circularly shifted coefficients.
}
\details{
  The "center of energy" technique of Wickerhauser and Hess-Nielsen
  (1996) is employed to find circular shifts for the wavelet
  sub-matrices such that the coefficients are aligned with the original
  series.  This corresponds to applying a (near) linear-phase filtering
  operation.
}
\references{
  Hess-Nielsen, N. and M. V. Wickerhauser (1996)
  Wavelets and time-frequency analysis,
  \emph{Proceedings of the IEEE}, \bold{84}, No. 4, 523-540.
  
  Percival, D. B. and A. T. Walden (2000)
  \emph{Wavelet Methods for Time Series Analysis},
  Cambridge University Press.
}
\seealso{
  \code{\link{phase.shift}}, \code{\link{modwt.2d}}.
}
\examples{
G1 <- G2 <- dnorm(seq(-512/4,512/4,length=512))
G <- 100 * zapsmall(G1 \%o\% G2)
J <- 6
G <- modwt.2d(G, "la8", J)
names.modwt <- c(t(sapply(c("LH","HL","HH"), paste, 1:J, sep="")), 
                 paste("LL", J, sep=""))
n <- 50
par(mfrow=c(3,3), mar=c(5,4,4,2)/1.9, pty="s")
for(i in names.modwt[1:9])
  image(G[[i]][256 + -n:n, 256 + -n:n], axes=FALSE, col=rainbow(128), main=i)
Gs <- shift.2d(G, "la8")
for(i in names.modwt[1:9])
  image(Gs[[i]][256 + -n:n, 256 + -n:n], axes=FALSE, col=rainbow(128), main=i)
}
\author{B. Whitcher}
\keyword{ts}
