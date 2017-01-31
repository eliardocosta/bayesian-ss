hpdPearsonVI <- function(len, kappa, psi, phi, w) {
  fun <- function(x) {
    (kappa - 1)*log(1 + len/x) - (kappa + psi)*log(1 + w*len/(phi + w*x))
  }
  sol <- uniroot(fun, lower = 1E-1, upper = 1E2, extendInt = "downX")
  return(sol$root)
}
