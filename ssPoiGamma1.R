ssPoiGamma1 <- function(crit, lam0, theta0, w, rho, len = NULL, 
                         len.max = NULL, R = 1E3, n0 = 0) {
  cat("\nCall for Poisson-gamma \n")
  cat("Criterion =", crit, "\n")
  if (crit == "CVM") cat("lam0 =", lam0, "; theta0 =", theta0, "; w =", w, "; eps =", eps, "\n")
  if (crit == "CCM1") cat("lam0 =", lam0, "; theta0 =", theta0, "; w =", w, "; l =", len, "\n")
  if (crit == "CCM2" || crit == "Aprox") cat("lam0 =", lam0, "; theta0 =", theta0, "; w =", w, "; l.max =", len.max, "\n")
  if (crit == "Aprox") {
    zrho <- qnorm(1 - rho/2)
    n <- (theta0/(w*lam0))*(((lam0/theta0)*(2*zrho/len.max)*(gamma(theta0 + 0.5)/gamma(theta0)))^2 - 1)
    cat("n =", ceiling(n), "\n")
  }
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
      }
    }
    cat("n =", n, "; cob. estimada =", mean(cov), "\n")
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
        len <- append(len, qgamma(1 - rho/2, shape = kappa,
                                       rate = psi) - qgamma(rho/2, 
                                                              shape = kappa, rate = psi))
      }
    }
    cat("n =", n, "; comp. estimado =", mean(len), "\n")
  } # FIM CRITERIO CCM2
} # FIM
