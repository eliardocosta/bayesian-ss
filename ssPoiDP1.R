ssPoiDP1 <- function(lam0, phi, w, rho, alpha, len.max, inc1 = 1E2 , inc2 = 5E1, R = 5E1) {
  cat("Call for PoiDP \n")
  cat("phi =", phi," ;", "w =", w, " ;", "alpha =", alpha ,"len.max =", len.max, "\n")
  n <- 0
  len.med <- len.max+1
  while (len.med > len.max) {
    n <- n+inc1
    #print(n)
    len <- numeric()
    for (i in 1:R) {
      obj.ppriori <- pprioriDP(n = n, alpha = alpha, lam0 = lam0, phi = phi)
      obj.dpriori <- extract.prob(obj.ppriori)
      #plot(obj.dpriori$lam, obj.dpriori$ddist, type = "h")
      lam <- sample(obj.dpriori$lam, n, replace = TRUE, prob = obj.dpriori$ddist)
      x <- rpois(n, w*lam)
      obj.ppost <- ppostPoiDP(x = x, w = w, phi = phi, lam0 = lam0, alpha = alpha, nburn = 5E2, nsam = 1E2)
      obj.dpost <- extract.prob(obj.ppost)
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1-rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
      }
      len <- append(len, length(conj.vals)*obj.ppriori$cgrid)
      cat(mean(len), " ")
    }
    cat("len.med e n \n")
    len.med <- mean(len)#; 
    cat(len.med, n, "\n")
  }
  # SEGUNDO LOOP COM INCREMENTO IGUAL A inc2
  n <- n-inc1
  len.med <- len.max+1
  while (len.med > len.max) {
    n <- n+inc2
    #print(n)
    len <- numeric()
    for (i in 1:R) {
      obj.ppriori <- pprioriDP(n = n, alpha = alpha, lam0 = lam0, phi = phi)
      obj.dpriori <- extract.prob(obj.ppriori)
      #plot(obj.dpriori$lam, obj.dpriori$ddist, type = "h")
      lam <- sample(obj.dpriori$lam, n, replace = TRUE, prob = obj.dpriori$ddist)
      x <- rpois(n, w*lam)
      obj.ppost <- ppostPoiDP(x = x, w = w, phi = phi, lam0 = lam0, alpha = alpha, nburn = 5E2, nsam = 1E2)
      obj.dpost <- extract.prob(obj.ppost)
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1-rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
      }
      len <- append(len, length(conj.vals)*obj.ppriori$cgrid)
      cat(mean(len), " ")
    }
    cat("len.med e n \n")
    len.med <- mean(len)#; 
    cat(len.med, n, "\n")
  }
  # TERCEIRO LOOP COM INCREMENTO IGUAL A 1
  n <- n-inc2
  len.med <- len.max+1
  while (len.med > len.max) {
    n <- n+1
    #print(n)
    len <- numeric()
    for (i in 1:R) {
      obj.ppriori <- pprioriDP(n = n, alpha = alpha, lam0 = lam0, phi = phi)
      obj.dpriori <- extract.prob(obj.ppriori)
      #plot(obj.dpriori$lam, obj.dpriori$ddist, type = "h")
      lam <- sample(obj.dpriori$lam, n, replace = TRUE, prob = obj.dpriori$ddist)
      x <- rpois(n, w*lam)
      obj.ppost <- ppostPoiDP(x = x, w = w, phi = phi, lam0 = lam0, alpha = alpha, nburn = 5E2, nsam = 1E2)
      obj.dpost <- extract.prob(obj.ppost)
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1-rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
      }
      len <- append(len, length(conj.vals)*obj.ppriori$cgrid)
      cat(mean(len), " ")
    }
    cat("len.med e n \n")
    len.med <- mean(len)#; 
    cat(len.med, n, "\n")
  }
  cat("n (PoiDP) =", n, "Comp est =", len.med, "\n")
} # FIM

# ssPoiDP1(lam0=10, phi=5, w=1, rho=0.05, alpha=0.5, len.max=2, R=100)
