ssPoiDP2 <- function(lam0, phi, w, rho, alpha, len.max, cgrid, inc1 = 1E2 , inc2 = 1E1, inc3 = 5, R = 1E2, n0 = 2) {
  cat("\nCall for PoiDP \n")
  cat("phi =", phi,"; w =", w, "; alpha =", alpha,"; len.max =", len.max, "\n")
  n <- n0
  len.med <- len.max + 1
  while (len.med > len.max) {
    n <- n + inc1
    len <- numeric()
    for (i in 1:R) {
      x <- rnbinom(n, mu = w*lam0, size = phi)
      obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
      obj.dpost <- extract.prob(obj.ppost)
  #   plot(obj.dpost$lam, obj.dpost$ddist, type = "h")
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist#/sum(obj.dpost$ddist)
      #soma1 <- sum(probs)
      #cat("soma igual a ", soma1, "\n")
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1 - rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
       # cat(conj.vals, "\n")
        #cat(probs, "\n")
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
          #cob[i] <- probs[vals == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
        #soma2 <- sum(cob)
        #cat(soma2, "\n")
      }
      len <- append(len, length(conj.vals)*obj.ppost$cgrid)
      #cat(mean(len), " ")
    }
    #cat("\nlen.med e n \n")
    len.med <- mean(len)#; 
    #cat(len.med, n, "\n")
  }
  # SEGUNDO LOOP COM INCREMENTO IGUAL A inc2
  n <- n - inc1
  len.med <- len.max + 1
  while (len.med > len.max) {
    n <- n + inc2
    len <- numeric()
    for (i in 1:R) {
      x <- rnbinom(n, mu = w*lam0, size = phi)
      obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
      obj.dpost <- extract.prob(obj.ppost)
   #   plot(obj.dpost$lam, obj.dpost$ddist, type = "h")
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist#/sum(obj.dpost$ddist)
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1 - rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
      }
      len <- append(len, length(conj.vals)*obj.ppost$cgrid)
      #cat(mean(len), " ")
    }
    #cat("\nlen.med e n \n")
    len.med <- mean(len)#; 
    #cat(len.med, n, "\n")
  }
  # TERCEIRO LOOP COM INCREMENTO IGUAL A inc3
  n <- n - inc2
  len.med <- len.max + 1
  while (len.med > len.max) {
    n <- n + inc3
    len <- numeric()
    for (i in 1:R) {
      x <- rnbinom(n, mu = w*lam0, size = phi)
      obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
      obj.dpost <- extract.prob(obj.ppost)
  #    plot(obj.dpost$lam, obj.dpost$ddist, type = "h")
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist#/sum(obj.dpost$ddist)
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1 - rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
      }
      len <- append(len, length(conj.vals)*obj.ppost$cgrid)
      #cat(mean(len), " ")
    }
    #cat("\nlen.med e n \n")
    len.med <- mean(len)#; 
    #cat(len.med, n, "\n")
  }
  # QUARTO LOOP COM INCREMENTO IGUAL A 1
  n <- n - inc3
  len.med <- len.max + 1
  while (len.med > len.max) {
    n <- n+1
    #print(n)
    len <- numeric()
    for (i in 1:R) {
      x <- rnbinom(n, mu = w*lam0, size = phi)
      obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
      obj.dpost <- extract.prob(obj.ppost)
   #   plot(obj.dpost$lam, obj.dpost$ddist, type = "h")
      vals <- obj.dpost$lam
      probs <- obj.dpost$ddist#/sum(obj.dpost$ddist)
      conj.vals <- numeric()
      cob <- 0
      while (sum(cob) < 1 - rho) {
        conj.ind <- which.max(probs) 
        conj.vals <- append(conj.vals, vals[conj.ind])
        for (i in 1:length(conj.vals)) {
          cob[i] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[i]]
        }
        vals <- vals[-conj.ind]
        probs <- probs[-conj.ind]
      }
      len <- append(len, length(conj.vals)*obj.ppost$cgrid)
      #cat(mean(len), " ")
    }
    #cat("\nlen.med e n \n")
    len.med <- mean(len)#; 
    #cat(len.med, n, "\n")
  }
  cat("n (PoiDP) =", n, "Comp est =", len.med, "\n")
}  # FIM
  
# ssPoiDP2(lam0 = lam0, phi = 5, w = w, rho = 0.05, alpha = 1, len.max = 4, inc1 = 1E2, inc2 = 5E1, inc3 = 1E1)