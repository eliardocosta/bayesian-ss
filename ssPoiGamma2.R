# SAMPLE SIZE USING POISSON-GAMMA MODEL UNDER DECISION 
# THEORETIC APPROACH
# LOSS FUNCTION 1: rho*tau + (a - lambda)^+ + (lambda-b)^+ -> lf = 1
# LOSS FUNCTION 2: gam*tau + (lam - m)^2/tau               -> lf = 2
ssPoiGamma2 <- function(lf, lam0, theta0, w, c, rho = NULL,
                         gam = NULL, nmax = 1E2, nrep = 1E1, lrep = 1E2,
                         plot = TRUE, ...) {
  cat("\nCall for Poisson-gamma \n")
  cat("Loss function =", lf, "\n")
  if (lf == 1) cat("lam0 =", lam0, "; theta0 =", theta0,"; w =", w, "; rho =", rho, "; c =", c, "\n")
  if (lf == 2) cat("lam0 =", lam0, "; theta0 =", theta0,"; w =", w, "; gam =", gam, "; c =", c, "\n")
  ns <- seq(1, nmax, by = 5)
  risk <- numeric()
  if (lf == 1) {
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        sn <- rnbinom(lrep, mu = n*w*lam0, size = n*theta0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        a <- qgamma(rho/2, shape = kappa, rate = psi)
        b <- qgamma(1 - rho/2, shape = kappa, rate = psi)
        loss <- (kappa/psi)*(pgamma(b, shape = kappa + 1, rate = psi, lower.tail = FALSE) - pgamma(a, shape = kappa + 1, rate = psi)) + c*n
        risk <- append(risk, mean(loss))
      }
    }
  } else if (lf == 2){
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        sn <- rnbinom(lrep, mu = n*w*lam0, size = n*theta0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        loss <- 2*sqrt(gam*kappa)/psi + c*n
        risk <- append(risk, mean(loss))
      }
    }
  }
  Y <- log(risk - c*rep(ns, each = nrep))
  mod <- lm(Y ~ I(log(rep(ns + 1, each = nrep))))
  E <- as.numeric(exp(mod$coef[1]))
  G <- as.numeric(-mod$coef[2])
  nmin <- ceiling((E*G/c)^(1/(G + 1))-1)
  if (plot == TRUE) {
    plot(rep(ns, each = nrep), risk, xlim = c(0, nmax), xlab = "n", ylab = "CT(n)")
    curve <- function(x) {c*x + E/(1 + x)^G}
    plot(function(x)curve(x), 0, nmax, col = "blue", add = TRUE)
    abline(v = nmin, col = "red")
  }
  cat("n =", nmin)
}
