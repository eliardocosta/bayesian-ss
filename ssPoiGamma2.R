# SAMPLE SIZE USING POISSON-GAMMA MODEL UNDER DECISION 
# THEORETIC APPROACH
# LOSS FUNCTION 1: rho*tau+(a-lambda)^+ + (lambda-b)^+ -> lf = 1
# LOSS FUNCTION 2: gam*tau+(lam-m)^2/tau               -> lf = 2
ssPoiGamma2 <- function(lf, lam0, theta0, w, c, rho = NULL,
                         gam = NULL, nmax = 1E3, nrep = 1E2, lrep = 1E2,
                         plot = TRUE, ...) {
  risk <- numeric()
  if (lf == 1) {
    for (n in 1:nmax) {
      for (i in 1:nrep) {
        loss <- numeric()
        #lam <- rgamma(lrep, shape = theta0, rate = theta0/lam0)
        #sn <- rpois(lrep, n*w*lam)
        sn <- rnbinom(lrep, mu = n*w*lam0, size = n*theta0)
        #sn <- rpois(lrep, n*w*lam0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        a <- qgamma(rho/2, shape = kappa, rate = psi)
        b <- qgamma(1 - rho/2, shape = kappa, rate = psi)
        tau <- (b - a)/2
        #loss <- rho*tau+ifelse(a-lam>0,a-lam,0)+ifelse(lam-b>0,lam-b,0)+c*n
        loss <- (kappa/psi)*(pgamma(b, shape = kappa + 1, rate = psi, lower.tail = FALSE) - pgamma(a, shape = kappa + 1, rate = psi)) + c*n
        risk <- append(risk, mean(loss))
      }
    }
  } else if (lf == 2){
    for (n in 1:nmax) {
      for (i in 1:nrep) {
        loss <- numeric()
        #lam <- rgamma(lrep, shape = phi, rate = phi/lam0)
        sn <- rnbinom(lrep, mu = n*w*lam0, size = n*theta0)
        #s <- rpois(lrep, n*w*lam)
        #for (j in 1:100) {
        #s <- rpois(1, n*w*lam0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        a <- kappa/psi - sqrt(kappa/gam)/psi
        b <- kappa/psi + sqrt(kappa/gam)/psi
        tau <- (b - a)/2
        m <- (a + b)/2
        #lam <- rgamma(lrep, shape = psi, rate = kappa)
        #loss <- append(loss, mean(gam*tau+(lam-m)^2/tau+c*n))
        loss <- 2*sqrt(gam*kappa)/psi + c*n
        #}
        risk <- append(risk, mean(loss))
      }
    }
  }
  Y <- log(risk - c*rep(1:nmax, each = nrep))
  mod <- lm(Y ~ I(log(rep(1:nmax + 1, each = nrep))))
  e <- as.numeric(exp(mod$coef[1]))
  g <- as.numeric(-mod$coef[2])
  nmin <- ceiling((e*g/c)^(1/(g + 1))-1)
  if (plot == TRUE) {
    plot(rep(1:nmax, each = nrep), risk, xlim = c(0, nmin + 0.5*nmax), xlab="n")
    curve <- function(x) {c*x + e/(1 + x)^g}
    plot(function(x)curve(x), 0, nmin + 0.5*nmax, col = "blue", add=TRUE)
    abline(v = nmin, col = "red")
  }
  cat("\nCall for Poisson-gamma \n")
  if (lf == 1) cat("lam0 =", lam0, ";theta0 =", theta0,"; w =", w, "; rho =", rho, "; c = ", c, "\n")
  if (lf == 2) cat("lam0 =", lam0, ";theta0 =", theta0,"; w =", w, "; gam =", gam, "; c = ", c, "\n")
  cat("n =", n)
}

# EXAMPLES
#
ssPoiGamma2(lf = 2, lam0=10, theta0 = 5, w = 0.5, c=0.01, gam=1/1, plot = FALSE, nmax = 50)
ssPoiGamma2(lf = 1, lam0=10, theta0 = 0.5, w = 1, c=0.01, rho = 0.1, plot = TRUE, nmax = 50)
