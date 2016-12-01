# Distribuição a priori usando a representação stick-breaking
#
pprioriDP <- function(n, alpha, lam0, phi, cgrid = 0.1) {
  #M <- ceiling(1-alpha*log(eps/(4*n)))
  M <- 1E2
  V <- rbeta(M-1, 1, alpha)
  V <- append(V, 1)
  p <- V[1]
  for (i in 2:M) {
    p[i] <- V[i]*(1-V[i-1])*p[i-1]/V[i-1]
  }
  lam <- rgamma(M, shape = phi, rate = phi/lam0)
  grid.val1 <- seq(0, ceiling(max(lam)), cgrid) # grid de valores que define a particao sobre possiveis valores de lambda
  #grid.val1 <- grid.val1[-1]
  grid.val1 <- tail(grid.val1, -1)
  ppriori <- numeric()
  #matrix(NA, nsam, length(grid.val)) # prob. acumulada a posteriori em relacao aos valores do gride por iteracao
  for (i in 1:length(grid.val1)) {
    ppriori[i] <- sum(p[which(lam <= grid.val1[i])])
  }
  #grid.val2 <- seq(0, ceiling(max(lam)), 0.5)
  #grid.val3 <- filter(grid.val2, rep(1/2, 2)) # pegando o ponto médio do intervalo
  #grid.val3 <- as.numeric(grid.val3[-length(grid.val3)])
  #dpriori <- numeric()
  #dpriori <- ppriori[1]
  #for (i in 2:length(ppriori)) {
  #  dpriori[i] <- ppriori[i]-ppriori[i-1]
  #}
  return(list(lam = grid.val1, fdist = ppriori, cgrid = cgrid))
} # FIM 
# pprioriDP(n = 5, alpha = 2, lam0 = 13, phi = 5)
