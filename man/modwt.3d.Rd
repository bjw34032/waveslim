% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/three_D.R
\name{modwt.3d}
\alias{modwt.3d}
\alias{imodwt.3d}
\title{Three Dimensional Separable Maximal Ovelrap Discrete Wavelet Transform}
\usage{
modwt.3d(x, wf, J = 4, boundary = "periodic")

imodwt.3d(y)
}
\arguments{
\item{x}{input array}

\item{wf}{name of the wavelet filter to use in the decomposition}

\item{J}{depth of the decomposition}

\item{boundary}{only \code{"periodic"} is currently implemented}

\item{y}{an object of class \code{modwt.3d}}
}
\description{
Three-dimensional separable maximal overlap discrete wavelet transform
(MODWT).
}
\author{
B. Whitcher
}
\keyword{ts}
