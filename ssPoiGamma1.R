ssPoiGamma1 <- function(lam0, theta0, w, rho, crit, len = NULL, 
                         len.max = NULL, R = 1E3, srep = 1E2, n0 = 0) {
  cat("\nCall for Poisson-gamma \n")
  if (crit == "CVM") cat("lam0 =", lam0,";theta0 =", theta0,"; w =", w, "; eps =", eps, "\n")
  if (crit == "CCM1") cat("lam0 =", lam0,";theta0 =", theta0,"; w =", w, "; l =", len, "\n")
  if (crit == "CCM2") cat("lam0 =", lam0,";theta0 =", theta0,"; w =", w, "; l.max =", len.max, "\n")
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
        sn <- rnbinom(1, mu = n*w*lam0, size = n*theta0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        a <- len/(exp((psi*len)/(kappa - 1)) - 1)
        cov <- append(cov, pgamma(a + len, shape = kappa, rate = psi) - pgamma(a, 
                                                                                    shape = kappa, rate = psi))
        #cov <- append(cov, mean(pgamma(a + len, shape = psi, rate = kappa) - pgamma(a, 
        #                                                                        shape = psi, rate = kappa)))
      }
    }
    cat("n =", n, "Cob est =", mean(cov), "\n")
  } # FIM CRITERIO CCM1
  if (crit == "CCM2") { # INICIO CRITERIO CCM2
    len <- len.max + 1
    n <- n0
    while (mean(len) > len.max) {
      n <- n + 1
      len <- numeric()
      for (i in 1:R) {
        sn <- rnbinom(1, mu = n*w*lam0, size = n*theta0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        #len <- append(len, mean(qgamma(1 - rho/2, shape = psi,
        #                               rate = kappa) - qgamma(rho/2, 
        #                                                    shape = psi, rate = kappa)))
        len <- append(len, qgamma(1 - rho/2, shape = kappa,
                                       rate = psi) - qgamma(rho/2, 
                                                              shape = kappa, rate = psi))
      }
    }
    cat("n =", n, "Comp est =", mean(len), "\n")
  } # FIM CRITERIO CCM2
} # FIM
