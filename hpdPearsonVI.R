hpdPearsonVI <- function(kappa, psi, phi, w, len = NULL, rho = NULL) {
  if (is.null(len)) {
    #require("rootSolve")
    #require("PearsonDS")
    fun <- function(x) c(F1 = ppearsonVI(x[2], a = kappa, b = psi, location = 0, scale = phi/w) - ppearsonVI(x[1], a = kappa, b = psi, location = 0, scale = phi/w) - 1 + rho,
                         F2 = dpearsonVI(x[1], a = kappa, b = psi, location = 0, scale = phi/w) - dpearsonVI(x[2], a = kappa, b = psi, location = 0, scale = phi/w))
    roots <- c(-2, -1)
    starts <- qpearsonVI(c(rho/2, 1 - rho/2), a = kappa, b = psi, location = 0, scale = phi/w)
    while (all(roots <= 0)) {
      sol.finder <- multiroot(f = fun, start = starts)
      roots <- sol.finder$root
      starts <- starts + 0.1
    }
  }
  if (is.null(rho)) {
    fun <- function(x) {
      (kappa - 1)*log(1 + len/x) - (kappa + psi)*log(1 + w*len/(phi + w*x))
    }
    sol <- uniroot(fun, lower = 1E-3, upper = 1E3, extendInt = "downX")
    roots <- c(sol$root, sol$root + len)
  }
  return(roots)
} # FIM
