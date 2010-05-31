\name{jumpsine}
\alias{jumpsine}
\title{Sine with Jumps at 0.625 and 0.875}
\description{
  \deqn{jumpsine(x) = 10\left( \sin(4{\pi}x) +
    I_{[0.625 < x \leq 0.875]}\right)}{%
    jumpsine(x) = 10*(sin(4*pi*x) + I_[0.625 < x <= 0.875])}
}
\usage{data(jumpsine)
}
\format{A vector containing 512 observations.
}
\source{S+WAVELETS.
}
\references{
  Bruce, A., and H.-Y. Gao (1996)
  \emph{Applied Wavelet Analysis with S-PLUS}, Springer: New York.
}
\keyword{datasets}
