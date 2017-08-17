#' @noRd
sumlist <- function(list){
   res <- 0
   for (i in seq(list)) res <- res + list[[i]]
   res
}

#' @noRd
par2Psi <- function(par, q){
   Psi <- {
      L <- diag(0, q)
      L[lower.tri(L, diag = TRUE)] <- par
      tcrossprod(L)
   }
   Psi
}

#' @noRd
glsfit <- function(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = TRUE){
   Sigmalist <- mapply(function(S, Z, na) 
      (S + Z %*% Psi %*% t(Z))[!na, !na, drop = FALSE], 
      Slist, Zlist, nalist, SIMPLIFY = FALSE)
   Ulist <- lapply(Sigmalist, chol)
   invUlist <- lapply(Ulist, function(U) backsolve(U, diag(ncol(U))))
   invtUXlist <- mapply(function(invU, X) crossprod(invU, X), 
                        invUlist, Xlist, SIMPLIFY = FALSE)
   invtUylist <- mapply(function(invU, y) crossprod(invU, y), 
                        invUlist, ylist, SIMPLIFY = FALSE)
   invtUX <- do.call("rbind", invtUXlist)
   invtUy <- do.call("rbind", invtUylist)
   coef <- as.numeric(qr.solve(invtUX, invtUy))
   if (onlycoef) 
      return(coef)
   list(coef = coef, Sigmalist = Sigmalist, Ulist = Ulist, invUlist = invUlist, 
        invtUXlist = invtUXlist, invtUylist = invtUylist, invtUX = invtUX, invtUy = invtUy)
}

#' @noRd
initpar <- function(Xlist, Zlist, ylist, Slist, nalist, q, control){
   initPsi <- diag(0.001, q)
   for (i in seq(control$igls.iter)) 
      initPsi <- iter.igls(initPsi, Xlist, Zlist, ylist, Slist, nalist, q)
   #  initPsi <- diag(0.001, p)
   return(vechMat(initPsi))
}

#' @noRd
gradchol.ml <- function(par, U, ind1, ind2, invSigmalist, reslist, Xlist, nalist, p){
   grad <- sapply(seq(length(par)), function(i) {
      A <- B <- C <- diag(0, p)
      A[ind2[i], ] <- B[, ind2[i]] <- U[ind1[i], ]
      C[ind2[i], ] <- C[, ind2[i]] <- 1
      D <- C * A + C * B
      gr <- sum(mapply(function(invSigma, res, na, X) {
         E <- crossprod(res, invSigma) %*% (X %*% D %*% t(X))[!na, !na] %*% 
            invSigma %*% res
         F <- sum(diag(invSigma %*% (X %*% D %*% t(X))[!na, !na]))
         return(as.numeric(0.5 * (E - F)))
      }, invSigmalist, reslist, nalist, Xlist))
   })
   grad
}

#' @noRd
gradchol.reml <- function(par, U, invtXWXtot, ind1, ind2, Xlist, invSigmalist, 
                          reslist, nalist, p){
   grad <- sapply(seq(length(par)), function(i) {
      A <- B <- C <- diag(0, p)
      A[ind2[i], ] <- B[, ind2[i]] <- U[ind1[i], ]
      C[ind2[i], ] <- C[, ind2[i]] <- 1
      D <- C * A + C * B
      gr <- sum(mapply(function(X, invSigma, res, na) {
         E <- invSigma %*% (X %*% D %*% t(X))[!na, !na] %*% invSigma
         F <- crossprod(res, E) %*% res
         G <- sum(diag(invSigma %*% (X %*% D %*% t(X))[!na, !na]))
         H <- sum(diag(invtXWXtot %*% crossprod(X, E) %*% 
                          X))
         return(as.numeric(0.5 * (F - G + H)))
      }, Xlist, invSigmalist, reslist, nalist))
   })
   grad
}


#' @noRd
# probably not need (cfr fracpol)
fracpol.eval <- function(x, p = c(1, 1), shift, scale, scaling = TRUE, xname = "x"){
   if (scaling){
      if (missing(shift)){
         shift <- 0
         if (length(x) > 1L){
            if (min(x) <= 0) {
               z <- diff(sort(x))
               shift <- min(z[z > 0]) - min(x)
               shift <- ceiling(shift * 10)/10
            }
         }
      }
      if (missing(scale)){
         scale <- 1
         if (length(x) > 1L){
            range <- mean(x + shift)
            scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
         }
      }
   } else {
      shift <- 0
      scale <- 1
   }
   scaling <- c(shift, scale)
   if (length(p) != 2L)
      stop("fracpol computes two-order fractional polynomials")
   x <- (x + scaling[1])/scaling[2]
   stopifnot(all(x > 0))
   
   Xp <- rbind(sapply(p, function(p) x^p))
   colnames(Xp) <- c(paste(xname, "^", p, sep = ""))
   if (any(p==0)){
      Xp[, p==0] <- log(x)
      colnames(Xp)[p==0] <- paste("log(", xname, ")", sep = "")
   }
   if (p[1] == p[2]){
      Xp[, 2] <- Xp[, 2]*log(x)
      colnames(Xp)[2] <- paste("log(", xname, ")", xname, "^", p[2], sep = "")
      if (any(p ==0))
         colnames(Xp)[2] <- paste("log(", xname, ")^2", sep = "")
   }
   attr(Xp, "p") <- p
   attr(Xp, "scaling") <- scaling
   Xp
}
