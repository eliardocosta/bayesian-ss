ssBNPearson1 <- function(lam0, theta0, phi, w, rho, crit, len = NULL, 
                        len.max = NULL, R1 = 1E2, R2 = 1E2, n0 = 0) {
  #install.packages("PearsonDS")
  #require("PearsonDS")
  #phi = 5
  #w = 1
  #theta0 = 5
  #lam0 = w*10/phi
  cat("\nCall for BNPearson \n")
  if (crit == "CVM") cat("phi =", phi,"; w =", w, "; theta0 =", theta0, "; eps =", eps, "\n")
  if (crit == "CCM1") cat("phi =", phi,"; w =", w, "; theta0 =", theta0, "; l =", len, "\n")
  if (crit == "CCM2") cat("phi =", phi,"; w =", w, "; theta0 =", theta0, "; l.max =", len.max, "\n")
  if (crit == "CVM") {
    break
  }
  if (crit == "CCM1") { # INICIO CRITERIO CCM1
    cov <- 0 
    n <- n0
    while (mean(cov) < 1 - rho) {
      n <- n + 1
      print(n)
      cov <- numeric()
      probs <- numeric()
      for (i in 1:R1) {
        lam <- rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, scale = phi/w)
        x <- rnbinom(length(lam), mu = w*w*lam/phi, size = phi)
        s <- sum(x)
        for (i in 1:R2) {
          #s <- rnbinom(srep, mu = n*w*lam0, size = n*phi)
          kappa <- theta0 + s
          psi <- theta0/lam0 + n*phi + 1
          a <- hpdPearsonVI(len = len, kappa = kappa, psi = psi, phi = phi, w = w)
          probs <- append(probs, 
                          ppearsonVI(a + len, a = kappa, b = psi, location = 0, 
                                     scale = phi/w) - ppearsonVI(a, a = kappa, b = psi, 
                                                             location = 0, scale = phi/w))
     #     print(probs)
        }
        cov <- append(cov, mean(probs))
      }
      print(mean(cov))
    }
    cat("n (BNPearson) =", n, "Cob est =", mean(cov), "\n")
  } # FIM CRITERIO CCM1
  if (crit == "CCM2") { # INICIO CRITERIO CCM2
    len <- len.max + 1
    n <- n0
    while (mean(len) > len.max) {
      n <- n + 1
      print(n)
      len <- numeric()
      lens <- numeric()
      for (i in 1:R1) {
        lam <- rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, scale = phi/w)
        x <- rnbinom(length(lam), mu = w*w*lam/phi, size = phi)
        s <- sum(x)
        for (i in 1:R2) {
          kappa <- theta0 + s
          psi <- theta0/lam0 + n*phi + 1
          lens <- append(lens, qpearsonVI(1 - rho/2, a = kappa, b = psi, location = 0, 
                                          scale = phi/w) - qpearsonVI(rho/2, a = kappa, 
                                                                      b = psi, location = 0,
                                                                      scale = phi/w))
        }
        len <- append(len, mean(lens))
      }
      print(mean(len))
    }
    cat("n (BNPearson) =", n, "Comp est =", mean(len), "\n")
  } # FIM CRITERIO CCM2
} # FIM
