ssPoiDP1 <- function(lam0, phi, w, rho, alpha, crit, len.max = NULL, len = NULL, eps = NULL,
                     cgrid = 1E-2, R = 1E2, n0 = 2, inc = c(1E2, 1E1, 5)) {
  cat("\nCall for PoiDP \n")
  if (crit == "CVM") cat("phi =", phi,"; w =", w, "; alpha =", alpha,"; eps =", eps, "\n")
  if (crit == "CCM1") cat("phi =", phi,"; w =", w, "; alpha =", alpha,"; l =", len, "\n")
  if (crit == "CCM2") cat("phi =", phi,"; w =", w, "; alpha =", alpha,"; l.max =", len.max, "\n")
  inc <- c(sort(inc, decreasing = TRUE), 1) 
  if (crit == "CVM") { # INICIO CRITERIO CVM
    n <- n0
    for (i in 1:length(inc)) {
      var.post <- eps + 1
      while (median(var.post) > eps) {
        n <- n + inc[i]
        var.post <- numeric()
        for (j in 1:R) {
          x <- rnbinom(n, mu = w*lam0, size = phi)
          obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
          obj.dpost <- extract.prob(obj.ppost)
          vals <- obj.dpost$lam
          probs <- obj.dpost$ddist
          med.post <- sum(vals*probs)
          var.post <- append(var.post, sum(probs*(vals - med.post)^2))
          cat(median(var.post), " ")
        }
        cat("\n cob.med e n \n")
        cat(median(var.post), n, "\n")
      }
      if (i < length(inc)) n <- n - inc[i]
    }
    cat("n (PoiDP) =", n, "Cob est =", median(var.post), "\n")
  } # FIM CRITERIO CVM
  if (crit == "CCM1") { # INICIO CRITERIO CCM1
    n <- n0
    for (i in 1:length(inc)) {
      cob <- 0
      while (mean(cob) < 1 - rho) {
        n <- n + inc[i]
        cob <- numeric()
        for (j in 1:R) {
          x <- rnbinom(n, mu = w*lam0, size = phi)
          obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
          obj.dpost <- extract.prob(obj.ppost)
          vals <- obj.dpost$lam
          probs <- obj.dpost$ddist
          conj.vals <- numeric()
          conj.ind <- numeric()
          for (k in 1:ceiling(len/obj.ppost$cgrid)) {
            ind.max <- which.max(probs) 
            conj.vals <- append(conj.vals, vals[ind.max])
            conj.ind <- append(conj.ind, which(obj.dpost$lam == vals[ind.max]))
            vals <- vals[-ind.max]
            probs <- probs[-ind.max]
          }
          cob <- append(cob, sum(obj.dpost$ddist[conj.ind]))
        }
      }
      if (i < length(inc)) n <- n - inc[i]
    }
    cat("n (PoiDP) =", n, "Cob est =", mean(cob), "\n")
  } # FIM CRITERIO CCM1
  if (crit == "CCM2") { # INICIO CRITÉRIO CCM2
    n <- n0
    for (i in 1:length(inc)) {
      len <- len.max + 1
      while (mean(len) > len.max) {
        n <- n + inc[i]
        len <- numeric()
        for (j in 1:R) {
          x <- rnbinom(n, mu = w*lam0, size = phi)
          obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, phi = phi, alpha = alpha, cgrid = cgrid, nsam = 5E1)
          obj.dpost <- extract.prob(obj.ppost)
          vals <- obj.dpost$lam
          probs <- obj.dpost$ddist
          conj.vals <- numeric()
          cob <- 0
          while (sum(cob) < 1 - rho) {
            ind.max <- which.max(probs) 
            conj.vals <- append(conj.vals, vals[ind.max])
            for (k in 1:length(conj.vals)) {
              cob[k] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[k]]
            }
            vals <- vals[-ind.max]
            probs <- probs[-ind.max]
          }
          len <- append(len, length(conj.vals)*obj.ppost$cgrid)
        }
      }
      if (i < length(inc)) n <- n - inc[i]
    }
    cat("n (PoiDP) =", n, "Comp est =", mean(len), "\n")
  } # FIM CRITERIO CCM2
}  # FIM
