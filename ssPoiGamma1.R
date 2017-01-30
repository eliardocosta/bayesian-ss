ssPoiGamma1 <- function(lam0, phi, w, rho, crit, len = NULL, 
                         len.max = NULL, R = 1E2, srep = 1E2, n0 = 0) {
  cat("\nCall for PoiGamma \n")
  if (crit == "CVM") cat("phi =", phi,"; w =", w, "; eps =", eps, "\n")
  if (crit == "CCM1") cat("phi =", phi,"; w =", w, "; l =", len, "\n")
  if (crit == "CCM2") cat("phi =", phi,"; w =", w, "; l.max =", len.max, "\n")
  if (crit == "CVM") {
    break
  }
  if (crit == "CCM1") { # INICIO CRITERIO CCM1
    cov <- 0 
    n <- n0
    while (mean(cov) < 1 - rho) {
      n <- n + 1
      cov <- numeric()
      for (i in 1:R) {
        s <- rnbinom(srep, mu = n*w*lam0, size = n*phi)
        psi <- phi + s
        kappa <- n*w + phi/lam0
        a <- len/(exp((kappa*len)/(psi - 1)) - 1)
        cov <- append(cov, mean(pgamma(a + len, shape = psi, rate = kappa) - pgamma(a, 
                                                                                shape = psi, rate = kappa)))
      }
    }
    cat("n (PoiGamma) =", n, "Cob est =", mean(cov), "\n")
  } # FIM CRITERIO CCM1
  if (crit == "CCM2") { # INICIO CRITERIO CCM2
    len <- len.max + 1
    n <- n0
    while (mean(len) > len.max) {
      n <- n + 1
      len <- numeric()
      for (i in 1:R) {
        s <- rnbinom(srep, mu = n*w*lam0, size = n*phi)
        psi <- phi + s
        kappa <- n*w + phi/lam0
        len <- append(len, mean(qgamma(1 - rho/2, shape = psi,
                                       rate = kappa) - qgamma(rho/2, 
                                                            shape = psi, rate = kappa)))
      }
    }
    cat("n (PoiGamma) =", n, "Comp est =", mean(len), "\n")
  } # FIM CRITERIO CCM2
} # FIM
