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
   if (!is.null(initPsi <- control$initPsi)) {
      if (is.vector(initPsi)) 
         initPsi <- xpndMat(initPsi)
      } else {
         initPsi <- diag(.0001, q)
   }
   if (control$igls.inititer != 0){
      for (i in seq(control$igls.inititer)){
         initPsi <- iter.igls(initPsi, Xlist, Zlist, ylist, Slist, nalist, q)
      } 
   }
   vechMat(t(chol(initPsi)))
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

#' @noRd
getGroups <- function (random, data) 
{
   if (is.null(random)) 
      return(matrix(seq(nrow(data))))
   random <- getList(random)
   groups <- lapply(random, function(form) {
      form[[2]] <- form[[2]][[3]]
      model.frame(form, data)[[1]]
   })
   groups[[1]] <- as.factor(groups[[1]])
   if ((len <- length(groups)) > 1L) 
      for (i in 2:len) groups[[i]] <- factor(paste(groups[[i - 
                                                              1]], groups[[i]], sep = "-"))
   groups <- do.call(cbind, lapply(groups, unclass))
   groups
}

#' @noRd
mixmeta.control <- function (optim = list(), showiter = FALSE, maxiter = 100, initPsi = NULL, 
          Psifix = NULL, Scor = NULL, addSlist = NULL, inputna = FALSE, 
          inputvar = 10^4, loglik.iter = "hybrid", igls.inititer = 10, 
          hessian = FALSE, vc.adj = TRUE, reltol = sqrt(.Machine$double.eps), 
          checkPD = NULL, set.negeigen = sqrt(.Machine$double.eps)) 
{
   optim <- modifyList(list(fnscale = -1, maxit = maxiter, reltol = reltol), 
                       optim)
   if (showiter) {
      optim$trace <- 6
      optim$REPORT <- 1
   }
   loglik.iter <- match.arg(loglik.iter, c("hybrid", "newton", 
                                           "igls", "rigls"))
   if (igls.inititer <= 0L) 
      igls.inititer <- 0
   list(optim = optim, showiter = showiter, maxiter = maxiter, 
        hessian = hessian, initPsi = initPsi, Psifix = Psifix, 
        Scor = Scor, addSlist = addSlist, inputna = inputna, 
        inputvar = inputvar, loglik.iter = loglik.iter, igls.inititer = igls.inititer, 
        vc.adj = vc.adj, reltol = reltol, checkPD = checkPD, 
        set.negeigen = set.negeigen)
}

#' @noRd
getList <- function (object) 
{
   if (is.list(object)) 
      object
   else list(object)
}
