# SAMPLE SIZE USING BN-PEARSON MODEL UNDER DECISION 
# THEORETIC APPROACH
# LOSS FUNCTION 1: rho*tau + (a - lambda)^+ + (lambda - b)^+ -> lf = 1
# LOSS FUNCTION 2: gam*tau + (lam - m)^2/tau                 -> lf = 2
ssBNPearson2 <- function(lf, lam0, theta0, phi, w, c, rho = NULL,
                        gam = NULL, nmax = 1E2, nrep = 1E1, lrep = 5E1,
                        plot = TRUE, ...) {
  cat("\nCall for BN-Pearson \n")
  cat("Loss function =", lf, "\n")
  if (lf == 1) cat("lam0 =", lam0, "; theta0 =", theta0, "; phi =", phi,"; w =", w, "; rho =", rho, "; c =", c, "\n")
  if (lf == 2) cat("lam0 =", lam0, "; theta0 =", theta0, "; phi =", phi,"; w =", w, "; gam =", gam, "; c =", c, "\n")
  ns <- seq(1, nmax, by = 5)
  risk <- numeric()
  if (lf == 1) {
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        sn <- numeric()
        for (i in 1:lrep) {
          lam <- rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, scale = phi/w)
          x <- rnbinom(length(lam), mu = w*lam, size = phi)
          sn <- append(sn, sum(x))
        }
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        a <- qpearsonVI(rho/2, a = kappa, b = psi, location = 0, scale = phi/w)
        b <- qpearsonVI(1 - rho/2, a = kappa, b = psi, location = 0, scale = phi/w)
        tau <- (b - a)/2
        loss <- (phi/w)*(ppearsonVI(b, a = kappa + 1, b = psi - 1, location = 0, scale = phi/w, lower.tail = FALSE) - ppearsonVI(a, a = kappa + 1, b = psi - 1, location = 0, scale = phi/w)) + c*n
        risk <- append(risk, mean(loss))
      }
    }
  } else if (lf == 2){
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        sn <- numeric()
        for (i in 1:lrep) {
          lam <- rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, scale = phi/w)
          x <- rnbinom(n, mu = w*lam, size = phi)
          sn <- append(sn, sum(x))
        }
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        qcon <- 1/(psi - 1)
        medpos <- (phi/w)*kappa/(psi - 1)
        varpos <- (kappa/(psi - 1)^2)*((medpos + 1)/(1 - qcon))*(phi/w)^2
        loss <- 2*sqrt(gam*varpos) + c*n
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
