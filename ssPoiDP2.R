# SAMPLE SIZE USING POISSON-DP MODEL UNDER DECISION 
# THEORETIC APPROACH (F0 = GAMMA)
# LOSS FUNCTION 3 
ssPoiDP2 <- function(lam0, theta0, alpha, w, c, nmax = 1E2, nrep = 1E1, plot = TRUE, ...) {
  cat("\nCall for Poisson-DP \n")
  cat("lam0 =", lam0, "; theta0 =", theta0, "; alpha =", alpha,"; w =", w, "; c =", c, "\n")
  ns <- seq(3, nmax, by = 5)
  risk <- numeric()
    for (n in ns) {  
      for (i in 1:nrep) {
        x <- rnbinom(n, mu = w*lam0, size = theta0)
        obj.vpost <- var.postDPmix(x = x, w = w, lam0 = lam0, theta0 = theta0, alpha = alpha)
        risk <- append(risk, sum(obj.vpost$varp*dnorm(obj.vpost$lam, mean = 10, sd = 1E1)) + c*n)
      }
    }
  Y <- log(risk - c*rep(ns, each = nrep))
  mod <- lm(Y ~ I(log(rep(ns + 1, each = nrep))))
  E <- as.numeric(exp(mod$coef[1]))
  G <- as.numeric(-mod$coef[2])
  nmin <- ceiling((E*G/c)^(1/(G + 1))-1)
  if (plot == TRUE) {
    plot(rep(ns, each = nrep), risk, xlim = c(0, nmin + 0.1*nmax), xlab = "n")
    curve <- function(x) {c*x + E/(1 + x)^G}
    plot(function(x)curve(x), 0, nmin + 0.1*nmax, col = "blue", add = TRUE)
    abline(v = nmin, col = "red")
  }
  cat("n =", nmin)
}
