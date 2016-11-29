# Distribuição a posteriori usando o algoritmo de ...
#
ppostPoiDP <- function (x, w, phi, lam0, alpha, nburn, nsam, cgrid = 0.1, eps = 1E-3) {
  nx <- length(x) # tamanho da amostra
  N <- min(nx, ceiling(1-alpha*log(eps/(4*nx)))) # num. de termos da representa??o em soma do DP
  lam <- rgamma(N, shape = phi, rate = phi/lam0) # valores iniciais para lambda
  p <- rep(1/N, N) # prob's iniciais p/ cada termo da soma truncada
  K <- numeric() # vetor indexador de cluster
  pstar <- matrix(NA, nx, N) 
  nclu <- numeric() # vetor com num. de cluster por itera??o
  mclu <- numeric() # vetor com m?dia do num. de cluster por itera??o
  for (r in 1:nburn) { # in?cio do loop do burn-in
    for(i in 1:nx) {
      pstar[i, ] <- p*dpois(x[i], w*lam)
    }
    for(i in 1:nx) {
      K[i] <- sample.int(N, size = 1, prob = pstar[i, ])
    }
    Ktable <- as.data.frame(table(K), stringsAsFactors = FALSE) # tabela freq. dos K's
    Kstar <- as.vector(Ktable[ ,1], mode="numeric") # vetor com os K's ?nicos
    nclu <- append(nclu, length(Kstar))
    mclu <- append(mclu, mean(nclu))
    m <- rep(0, N)
    for (j in 1:length(Ktable[ ,2])) {
      for (i in 1:N) {
        if (i == Ktable[j, 1]) m[i] <- Ktable[j, 2]
      }
    }
    D <- rgamma(length(m), shape = alpha/N + m, rate = 1)
    p <- D/sum(D) # atualizando vetor 'p'
    for (j in 1:length(Kstar)) {
      set = numeric()
      for (i in 1:nx) {
        if (K[i] == Kstar[j]) set = append(set, i)
      }
      lam[Kstar[j]] <- rgamma(1, shape = phi+sum(x[set]), 
                              rate = length(set)*w+phi/lam0)
    }
  } # fim do loop do burn-in
  # in?cio da amostragem da posteroiri
  pesos <- matrix(NA, nsam, N) # pesos a posteriori da soma truncada amostrados por itera??o
  lam.post <- matrix(NA, nsam, N) # lambda's a posteroiri amostrados por itera??o
  for (r in 1:nsam) { # in?cio do loop p/ amostras a posteriori
    for(i in 1:nx) {
      pstar[i, ] <- p*dpois(x[i], w*lam)
    }
    for(i in 1:nx) {
      K[i] <- sample.int(N, size = 1, prob = pstar[i, ])
    }
    Ktable <- as.data.frame(table(K), stringsAsFactors = FALSE) # tabela freq. dos K's
    Kstar <- as.vector(Ktable[ ,1], mode="numeric") # vetor com os K's ?nicos
    m <- rep(0, N)
    for (j in 1:length(Ktable[ ,2])) {
      for (i in 1:N) {
        if (i == Ktable[j, 1]) m[i] <- Ktable[j, 2]
      }
    }
    D <- rgamma(length(m), shape = alpha/N + m, rate = 1)
    p <- D/sum(D) # atualizando vetor 'p'
    for (j in 1:length(Kstar)) {
      set = numeric()
      for (i in 1:nx) {
        if (K[i] == Kstar[j]) set = append(set, i)
      }
      lam[Kstar[j]] <- rgamma(1, shape = phi+sum(x[set]), 
                              rate = length(set)*w+phi/lam0)
    }
    lam.post[r, ] <- lam
    pesos[r, ] <- p
  } # fim do loop p/ amostras a posteriori
  grid.val <- seq(0, ceiling(max(lam.post)), cgrid) # grid de valores que define a parti??o sobre poss?veis valores de lambda
  #grid.val <- grid.val[-1]
  grid.val <- tail(grid.val, -1)
  ppost.it <- matrix(NA, nsam, length(grid.val)) # prob. acumulada a posteriori em rela??o aos valores do gride por itera??o
  for (i in 1:nsam) {
    for (j in 1:length(grid.val)) {
      ppost.it[i, j] <- sum(p[which(lam.post[i, ] <= grid.val[j])])  
    }
  }
  ppost <- apply(ppost.it, 2, mean) # estimativa da prob. acumulada a posterirori
  return(list(lam = grid.val, fdist = ppost, mclu = mclu))
} # FIM

# ppostPoiDP(x = x, w = w, phi = phi, lam0 = lam0, alpha = alpha, nburn = 500, nsam = 1E3)