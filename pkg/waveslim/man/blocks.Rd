\name{blocks}
\alias{blocks}
\title{A Piecewise-Constant Function}
\description{
    \deqn{blocks(x) = \sum_{j=1}^{11}(1 + {\rm sign}(x-p_j)) h_j / 2}{%
      blocks(x) = sum[j=1,11] (1 + sign(x - p_j)) h_j/2}
}
\usage{data(blocks)
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
