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
        invtUXlist = invtUXlist, invtUX = invtUX, invtUy = invtUy)
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

#' Wald Test for Model Coefficients
#' 
#' @description Computes a Wald chi-squared test for 1 or more coefficients, given their variance-covariance matrix.
#' 
#' @param Sigma a var-cov matrix, usually extracted from one of the fitting functions.
#' @param b a vector of coefficients with var-cov matrix \code{Sigma}. These coefficients are usually extracted 
#' from one of the fitting functions available in \code{R}.
#' @param Terms an optional integer vector specifying which coefficients should be jointly tested, using a Wald 
#' chi-squared or F test. Its elements correspond to the columns or rows of the var-cov matrix given in \code{Sigma}. 
#' Default is \code{NULL}.
#' @param L an optional matrix conformable to \code{b}, such as its product with \code{b} gives 
#' the linear combinations of the coefficients to be tested. Default is \code{NULL}.
#' @param H0 a numeric vector giving the null hypothesis for the test. It must be as long as \code{Terms} or must 
#' have the same number of columns as \code{L}. Default to 0 for all the coefficients to be tested.
#' @param x Object of class "waldtest".
#' @param digits number of decimal places for displaying test results. Default to 2.
#' @param \dots further arguments passed to or from other methods.
#'
#' @details The \code{waldtest} and the method \code{print.waldtest} are taken from the \code{aod} package and 
#' simplified for ease of use.
#' 
#' @return An object of class \code{waldtest}, printed with \code{print.waldtest}.
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{aod}, \code{\link{summary.dosresmeta}}
#' 
#' @examples
#' ## Load data and run the model
#' data("alcohol_cvd")
#' model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                     se = se, cases = cases, n = n, data = alcohol_cvd) 
#'                     
#' ## Test for significance of the overall dose-response association
#' waldtest(b = coef(model), Sigma = vcov(model), Terms = 1:nrow(vcov(model)))
#'
#' @rdname waldtest
#' @export waldtest

waldtest <- function(Sigma, b, Terms = NULL, L = NULL, H0 = NULL){
   if (is.null(Terms) & is.null(L)) 
      stop("One of the arguments Terms or L must be used.")
   if (!is.null(Terms) & !is.null(L)) 
      stop("Only one of the arguments Terms or L must be used.")
   if (is.null(Terms)) {
      w <- nrow(L)
      Terms <- seq(length(b))[colSums(L) > 0]
   }
   else w <- length(Terms)
   if (is.null(H0)) 
      H0 <- rep(0, w)
   if (w != length(H0)) 
      stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
   if (is.null(L)) {
      L <- matrix(rep(0, length(b) * w), ncol = length(b))
      for (i in 1:w) L[i, Terms[i]] <- 1
   }
   dimnames(L) <- list(paste("L", as.character(seq(NROW(L))), 
                             sep = ""), names(b))
   f <- L %*% b
   V <- Sigma
   mat <- solve(L %*% V %*% t(L))
   stat <- t(f - H0) %*% mat %*% (f - H0)
   p <- 1 - pchisq(stat, df = w)
   res <- list(chitest = c(chi2 = stat, df = w, P = p), 
               Terms = Terms, b = b, H0 = H0)
   class(res) <- "waldtest"
   res   
}


#' @rdname waldtest
#' @method print waldtest
#' @export
print.waldtest <- function (x, digits = 2, ...){
   Terms <- x[["Terms"]]
   b <- x[["b"]]
   H0 <- x[["H0"]]
   v <- x[["chitest"]]
   namb <- names(b)[Terms]
   cat("Wald test:\n", "----------\n", sep = "")
   cat("\nChi-squared test:\n")
   cat("X2 = ", format(v["chi2"], digits = digits, nsmall = 1), 
       ", df = ", v["df"], ", P(> X2) = ", format(v["P"], digits = digits, 
                                                  nsmall = 1), "\n", sep = "")
}

#' Fractional Polynomials
#' 
#' @description Two-order fractional polynomials transformation for continuous covariates.
#' 
#' @param x a numeric vector.
#' @param p a vector of length 2 with the powers of x to be included.
#' @param scaling a logical indicating if the measurements are scaled prior to model fitting.
#' @param shift optional scalar representing the shift, if \code{scaling = TRUE}. If not specified it is se
#' internally equal to 0. 
#' @param scale optional scalar representing the scale, if \code{scaling = TRUE}. If not specified it is se
#' internally equal to 1.
#' 
#' @details The \code{fracpol} is based on the \code{FP} function in the \code{mboost} package.
#' See \code{help(FP)} for more details.
#' 
#' @return A matrix including the trasformations corresponding to the input values.
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{mboost}, \code{rcs.eval}
#' @references
#' 
#' Royston, Patrick, and Douglas G. Altman. "Regression using fractional polynomials of 
#' continuous covariates: parsimonious parametric modelling." Applied Statistics (1994): 429-467.
#' 
#' @examples
#' ## Load data and run the model
#' data("alcohol_cvd")
#' 
#' with(alcohol_cvd, fracpol(dose, p = c(.5, .5)))
#' 
#' model <- dosresmeta(formula = logrr ~ fracpol(dose, p = c(.5, .5)), type = type, id = id,
#'                     se = se, cases = cases, n = n, data = alcohol_cvd) 
#'                     
#' ## Test for significance of the overall dose-response association
#' waldtest(b = coef(model), Sigma = vcov(model), Terms = 1:nrow(vcov(model)))
#'
#' @export fracpol

fracpol <- function(x, p = c(1, 1), shift, scale, scaling = TRUE){
   xname <- deparse(substitute(x))

   if (scaling){
      if (missing(shift)){
         shift <- 0
         if (min(x) <= 0) {
            z <- diff(sort(x))
            shift <- min(z[z > 0]) - min(x)
            shift <- ceiling(shift * 10)/10
         }
      }
      if (missing(scale)){
         range <- mean(x + shift)
         scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
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
   
   Xp <- sapply(p, function(p) x^p)
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
