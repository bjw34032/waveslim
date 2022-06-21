#' Calculating Thomson's Spectral Multitapers by Inverse Iteration
#' 
#' The following function links the subroutines in "bell-p-w.o" to an R
#' function in order to compute discrete prolate spheroidal sequences (dpss).
#' 
#' Spectral estimation using a set of orthogonal tapers is becoming widely used
#' and appreciated in scientific research.  It produces direct spectral
#' estimates with more than 2 df at each Fourier frequency, resulting in
#' spectral estimators with reduced variance.  Computation of the orthogonal
#' tapers from the basic defining equation is difficult, however, due to the
#' instability of the calculations -- the eigenproblem is very poorly
#' conditioned.  In this article the severe numerical instability problems are
#' illustrated and then a technique for stable calculation of the tapers --
#' namely, inverse iteration -- is described. Each iteration involves the
#' solution of a matrix equation.  Because the matrix has Toeplitz form, the
#' Levinson recursions are used to rapidly solve the matrix equation.  FORTRAN
#' code for this method is available through the Statlib archive.  An
#' alternative stable method is also briefly reviewed.
#' 
#' @param n length of data taper(s)
#' @param k number of data tapers; 1, 2, 3, ... (do not use 0!)
#' @param nw product of length and half-bandwidth parameter (w)
#' @param nmax maximum possible taper length, necessary for FORTRAN code
#' @return \item{v}{matrix of data tapers (cols = tapers)}
#' \item{eigen}{eigenvalue associated with each data taper} \item{iter}{total
#' number of iterations performed} \item{n}{same as input}
#' \item{w}{half-bandwidth parameter} \item{ifault}{0 indicates success, see
#' documentation for "bell-p-w" for information on non-zero values}
#' @author B. Whitcher
#' @seealso \code{\link{sine.taper}}.
#' @references B. Bell, D. B. Percival, and A. T. Walden (1993) Calculating
#' Thomson's spectral multitapers by inverse iteration, \emph{Journal of
#' Computational and Graphical Statistics}, \bold{2}, No. 1, 119-130.
#' 
#' Percival, D. B. and A. T. Walden (1993) \emph{Spectral Estimation for
#' Physical Applications: Multitaper and Conventional Univariate Techniques},
#' Cambridge University Press.
#' @keywords ts
#' @export dpss.taper
dpss.taper <- function(n, k, nw = 4, nmax = 2^(ceiling(log(n,2)))) {
  if(n > nmax)
    stop("length of taper is greater than nmax")
  w <- nw/n
  if(w > 0.5)
    stop("half-bandwidth parameter (w) is greater than 1/2")
  if(k <= 0)
    stop("positive dpss order (k) required")
  v <- matrix(0, nrow = nmax, ncol = (k + 1))
  storage.mode(v) <- "double"
  out <- .Fortran(C_dpss,
                  nmax = as.integer(nmax),
                  kmax = as.integer(k),
                  n = as.integer(n),
                  w = as.double(w),
                  v = v,
                  sig = double(k + 1),
                  totit = integer(1),
                  sines = double(n),
                  vold = double(n),
                  u = double(n),
                  scr1 = double(n),
                  ifault = integer(1))
  ##list(v = out$v[1:n, 1:k], eigen = out$sig[-1] + 1, iter = 
  ##     out$totiTRUE, n = out$n, w = out$w, ifault = out$ifault)
  return(out$v[1:n, 1:k])
}



#' Computing Sinusoidal Data Tapers
#' 
#' Computes sinusoidal data tapers directly from equations.
#' 
#' See reference.
#' 
#' @param n length of data taper(s)
#' @param k number of data tapers
#' @return A vector or matrix of data tapers (cols = tapers).
#' @author B. Whitcher
#' @seealso \code{\link{dpss.taper}}.
#' @references Riedel, K. S. and A. Sidorenko (1995) Minimum bias multiple
#' taper spectral estimation, \emph{IEEE Transactions on Signal Processing},
#' \bold{43}, 188-195.
#' @keywords ts
#' @export sine.taper
sine.taper <- function(n, k) {
  tapers <- NULL
  for(i in 1:k)
    tapers <- cbind(tapers, sqrt(2/(n+1)) * sin((pi*i*1:n)/(n+1)))
  return(tapers)
}
