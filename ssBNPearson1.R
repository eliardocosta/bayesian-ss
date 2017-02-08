ssBNPearson1 <- function(crit, lam0, theta0, phi, w, rho, len = NULL, 
                        len.max = NULL, R = 1E3, n0 = 0) {
  cat("\nCall for BN-Pearson \n")
  cat("Criterion =", crit, "\n")
  if (crit == "CVM") cat("phi =", phi,"; w =", w, "; theta0 =", theta0, "; eps =", eps, "\n")
  if (crit == "CCM1") cat("lam0 =", lam0, "; theta0 =", theta0, "; phi =", phi, "; w =", w, "; l =", len, "\n")
  if (crit == "CCM2") cat("lam0 =", lam0, "; theta0 =", theta0, "; phi =", phi, "; w =", w, "; l.max =", len.max, "\n")
  if (crit == "CVM") {
    break
  }
  if (crit == "CCM1") { # INICIO CRITERIO CCM1
    cov <- 0 
    n <- n0
    while (mean(cov) < 1 - rho) {
      n <- n + 1
      cov <- numeric()
      probs <- numeric()
      for (i in 1:R) {
        lam <- rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, scale = phi/w)
        x <- rnbinom(n, mu = w*lam, size = phi)
        sn <- sum(x)
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        a <- hpdPearsonVI(len = len, kappa = kappa, psi = psi, phi = phi, w = w)
        cov <- append(cov, 
                          ppearsonVI(a + len, a = kappa, b = psi, location = 0, 
                                     scale = phi/w) - ppearsonVI(a, a = kappa, b = psi, 
                                                             location = 0, scale = phi/w))
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
        lam <- rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, scale = phi/w)
        x <- rnbinom(n, mu = w*lam, size = phi)
        sn <- sum(x)
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        len <- append(len, qpearsonVI(1 - rho/2, a = kappa, b = psi, location = 0, 
                                          scale = phi/w) - qpearsonVI(rho/2, a = kappa, 
                                                                      b = psi, location = 0,
                                                                      scale = phi/w))
      }
    }
    cat("n =", n, "; comp. estimado =", mean(len), "\n")
  } # FIM CRITERIO CCM2
} # FIM


ssBNPearson1(crit ="CCM1", lam0=10, theta0=11, phi=5, w=1, rho=0.05, len = 2)
