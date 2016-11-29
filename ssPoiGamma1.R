ssPoiGamma1 <- function(lam0, phi, w, rho, crit, len = NULL, 
                         len.max = NULL, R = 1E3, srep = 1E2) {
  cat("Call for PoiGamma \n")
  cat("phi =", phi," ;", "w =", w, " ;", "len.max =", len.max, "\n")
  if (crit == "ACC") {
    cov <- 0
    n <- 0
    while (mean(cov) < 1-rho) {
      #print(mean(cov))
      n <- n+1
      cov <- numeric()
      for (i in 1:R) {
        s <- rnbinom(srep, mu = n*w*lam0, size = n*phi)
        psi <- phi+s
        kappa <- n*w+phi/lam0
        a <- len/(exp((kappa*len)/(psi-1))-1)
        cov <- append(cov, mean(pgamma(a+len, shape = psi, rate = kappa)-pgamma(a, 
                                                                                shape = psi, rate = kappa)))
      }
    }
  }
  if (crit == "ALC") {
    len <- len.max+1
    n <- 0
    while (mean(len) > len.max) {
      #print(mean(len))
      n <- n+1
      len <- numeric()
      for (i in 1:R) {
        s <- rnbinom(srep, mu = n*w*lam0, size = n*phi)
        psi <- phi+s
        kappa <- n*w+phi/lam0
        len <- append(len, mean(qgamma(1-rho/2, shape = psi,
                                       rate = kappa)-qgamma(rho/2, 
                                                            shape = psi, rate = kappa)))
      }
      #print(mean(len))
    }
    cat("n (PoiGamma) =", n, "Comp est =", mean(len), "\n")
    cat("\n")
  }
  
  #return(n)
} # FIM