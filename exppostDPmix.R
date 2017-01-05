exp.postDPmix <- function(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = 1E-2, nsam = 1E2) {
  n <- length(x)
  samcon.lam <- rdistnu(nsam = nsam, x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha)
  grid <- seq(0, ceiling(max(samcon.lam)), cgrid)
  probs <- matrix(NA, nsam, length(grid))
  probs[ ,1] <- rep(0, nsam)
  for (i in 1:nsam) {
    for (j in 2:length(grid)) {
      probs[i, j] <- (alpha/(alpha + n))*pgamma(grid[j], shape = phi, rate = phi/lam0)+
        (1/(alpha + n))*(length(which(samcon.lam[i, ] <= grid[j])))
    }
  }
  return(list(lam = grid, fdist = apply(probs, 2, mean), cgrid = cgrid))
} # FIM
