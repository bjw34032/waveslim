hosking.sim <- function(n, acvs) {
  .C(C_hosking, tseries=rnorm(n), as.integer(n), as.double(acvs[1:n]))$tseries
}
