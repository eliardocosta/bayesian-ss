extract.prob <- function(obj) { # obter a distribuição de probabilidade do 'obj'
  grid.val <- filter(obj$lam, rep(1/2, 2)) # pegando o ponto médio do intervalo
  #grid.val <- as.numeric(grid.val[-length(grid.val)])
  grid.val <- as.numeric(head(grid.val, -1))
  dprob <- numeric()
  for (i in 2:length(obj$fdist)) {
    dprob[i] <- obj$fdist[i]-obj$fdist[i-1]
  }
  dprob <- dprob[-1]
  return(list(lam = grid.val, ddist = dprob))
} # FIM