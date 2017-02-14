var.postDPmix <- function(x = x, w = w, lam0 = lam0, theta0 = theta0, alpha = alpha, cgrid = 5E-2, nsam = 5E1) {
  n <- length(x)
  samcon.lam <- rdistnu(nsam = nsam, x = x, w = w, lam0 = lam0, theta0 = theta0, alpha = alpha)
  grid <- seq(0, ceiling(max(samcon.lam)), cgrid)
  vars <- matrix(NA, nsam, length(grid))
  vars[ ,1] <- rep(0, nsam)
  for (i in 1:nsam) {
    for (j in 2:length(grid)) {
      vars[i, j] <- ((alpha + n + 1)*(alpha + n)^2)^(-1)*(
        alpha^2*pgamma(grid[j], shape = theta0, rate = theta0/lam0)*pgamma(grid[j], shape = theta0, rate = theta0/lam0, lower.tail = FALSE) + 
        alpha*n*pgamma(grid[j], shape = theta0, rate = theta0/lam0) +
        length(which(samcon.lam[i, ] <= grid[j]))*(n - length(which(samcon.lam[i, ] <= grid[j]))) +
        alpha*(1 - 2*pgamma(grid[j], shape = theta0, rate = theta0/lam0))*length(which(samcon.lam[i, ] <= grid[j]))  
      )
    }
  }
  return(list(lam = grid, varp = apply(vars, 2, mean), cgrid = cgrid))
} # FIM 
