rdistnu <- function(nsam, x, w, lam0, phi, alpha, nburn = 1E2) {
  n <- length(x)
  lam <- rgamma(n, shape = phi, rate = phi/lam0)
    for (r in 1:nburn) { # COMECO BURN-IN
      for (i in 1:n) { 
      q0 <- alpha*dnbinom(x[i], mu = w*lam0, size = phi)
      qk <- dpois(x[i], lambda = w*lam[-i])
      cn <- q0 + sum(qk) # constante de normalizacao
      q0n <- q0/cn
      qkn <- qk/cn
      u <- runif(1)
      lam[i] <- ifelse(u <= q0n, rgamma(1, shape = phi + x[i], rate = (w + phi/lam0)),
                       sample(x = lam[-i], size = 1, prob = qkn))
      }
    } # FIM BURN-IN
    sam.lam <- matrix(NA, nsam, n)
    for (r in 1:nsam) { # COMECO AMOSTRAGEM
      for (i in 1:n) {
        q0 <- alpha*dnbinom(x[i], mu = w*lam0, size = phi)
        qk <- dpois(x[i], lambda = w*lam[-i])
        cn <- q0 + sum(qk) # constante de normalizacao
        q0n <- q0/cn
        qkn <- qk/cn
        u <- runif(1)
        lam[i] <- ifelse(u <= q0n, rgamma(1, shape = phi + x[i], rate = (w + phi/lam0)),
                         sample(x = lam[-i], size = 1, prob = qkn))
      }
      sam.lam[r, ] <- lam
    } # FIM AMOSTRAGEM
    return(sam.lam)
} # FIM
