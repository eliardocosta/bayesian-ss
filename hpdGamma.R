hpdGamma <- function(kappa, psi, rho = NULL, len = NULL) {
  if (is.null(len)) {
    require("rootSolve")
    fun <- function(x) c(F1 = pgamma(x[2], shape = kappa, rate = psi) - pgamma(x[1], shape = kappa, rate = psi) - 1 + rho,
                         F2 = dgamma(x[1], shape = kappa, rate = psi) - dgamma(x[2], shape = kappa, rate = psi))
    roots <- c(-2, -1)
    starts <- qgamma(c(rho/2, 1 - rho/2), shape = kappa, rate = psi)
    while (all(roots <= 0)) {
      sol.finder <- multiroot(f = fun, start = starts)
      roots <- sol.finder$root
      starts <- starts + 0.1
    }
  }
  if (is.null(rho)) {
    a <- len/(exp((psi*len)/(kappa - 1)) - 1)
    roots <- c(a, a + len)
  }
  return(roots)  
} # FIM
