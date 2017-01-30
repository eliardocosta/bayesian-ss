hpdPearsonVI <- function(len, kappa, psi, phi, w) {
  fun <- function(x) {
    (kappa - 1)*log(1 + len/x) - (kappa + psi)*log(1 + w*len/(phi + w*x))
  }
  sol <- uniroot(fun, lower = 0, upper = 1E5)
  return(sol$root)
}
