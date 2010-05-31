shift.2d <- function(z, wf, inv=FALSE)
{
  coe <- function(g) sum(0:(length(g)-1) * g^2) / sum(g^2)
  
  h <- wave.filter(wf)$hpf
  g <- wave.filter(wf)$lpf
  
  J <- (length(z) - 1) / 3
  n <- dim(z[[1]])[1]
  
  nu.H <- round(2^(1:J-1) * (coe(g) + coe(h)) - coe(g), 0)
  nu.H <- ifelse(nu.H/n < 1, nu.H, nu.H - trunc(nu.H/n) * n)
  nu.G <- round((2^(1:J) - 1) * coe(g), 0)
  nu.G <- ifelse(nu.G/n < 1, nu.G, nu.G - trunc(nu.G/n) * n)
  
  if(!inv) {
    ## Apply the phase shifts
    for(j in 0:(J-1)) {
      H.order <- c((nu.H[j+1]+1):n, 1:nu.H[j+1])
      G.order <- c((nu.G[j+1]+1):n, 1:nu.G[j+1])
      z[[3*j+1]] <- z[[3*j+1]][G.order, H.order]
      z[[3*j+2]] <- z[[3*j+2]][H.order, G.order]
      z[[3*j+3]] <- z[[3*j+3]][H.order, H.order]
    } 
    z[[3*J+1]] <- z[[3*J+1]][G.order, G.order]
  }
  else {
    ## Apply the phase shifts "reversed"
    for(j in 0:(J-1)) {
      H.order <- c((n-nu.H[j+1]+1):n, 1:(n-nu.H[j+1]))
      G.order <- c((n-nu.G[j+1]+1):n, 1:(n-nu.G[j+1]))
      z[[3*j+1]] <- z[[3*j+1]][G.order, H.order]
      z[[3*j+2]] <- z[[3*j+2]][H.order, G.order]
      z[[3*j+3]] <- z[[3*j+3]][H.order, H.order]
    }
    z[[3*J+1]] <- z[[3*J+1]][G.order, G.order]
  }
  return(z)
}

