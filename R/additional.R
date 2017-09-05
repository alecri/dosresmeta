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


#' Computes statistics to evaluate the goodness-of-fit from dosresmeta Objects
#' 
#' @name gof
#' 
#' @description This function computes statistics to evaluate the goodness-of-fit for dose-response meta-analysis.
#' It implements the deviance test, the coefficient of determination, and a dataframe useful for a decorrelated residuals-versus-exposure plot.
#' See reference for more details
#' 
#' @param object an object of class \code{dosresmeta} produced by \code{\link{dosresmeta}}.
#' @param x an object of class \code{gof.dosresmeta} produced by \code{gof}.
#' @param fixed logical for selecting fixed model. By default equal to \code{TRUE}.
#' @param digits an integer specifying the number of digits to which printed results must be rounded.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @return A list of class \code{gof.dosresmeta} containing the following
#' \tabular{ll}{
#' \code{tdata} \tab a dataframe with the decorrelated variables (y*, X*, e*).\cr
#' \code{R2} \tab Coefficient of determination R^2.\cr
#' \code{deviance} \tab Deviance test.\cr
#' }
#' 
#'
#' @examples
#' ## Loading the data
#' data("milk_ov")
#' 
#' ## Linear dose-response model
#' lin <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                  se = se, cases = case, n = n, data = milk_ov)
#'                  
#' ## Display goodness-of-fit statistics
#' gof(lin)
#' 
#' ## Meta-regression model
#' lin_reg <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'   se = se, cases = case, n = n, data = milk_ov,
#'   mod = ~ type)
#'
#' ## Display goodness-of-fit statistics for meta-regression model
#' gof(lin_reg)
#'  
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references 
#' 
#' Discacciati A, Crippa A, Orsini N. Goodness of fit tools for dose-response meta-analysis of binary outcomes. 
#' Research synthesis methods. 2015 Jan 1.
#' 
#' @export gof

gof <- function(object, fixed = TRUE){
   mf <- model.frame(object)
   mfmX <- object$model
   y <- as.matrix(model.response(mfmX, "numeric"))
   id <- object$id
   nay <- is.na(y)
   X <- model.matrix(attr(mfmX, "terms"), data = mfmX)
   X <- X[, -grep("(Intercept)", colnames(X)), drop = FALSE]
   Slist <- object$Slist
   v <- object$v
   dim <- object$dim
   vlist <- lapply(unique(id), function(j) cbind(v[id == j]))
   ylist <- lapply(unique(id), function(j) 
      y[id == j, , drop = FALSE][!nay[j, ], , drop = FALSE])
   Xlist <- lapply(unique(id),
                   function(j) 
                      X[id == j, , drop = FALSE][!nay[j, ], , drop = FALSE])
   nalist <- lapply(unique(id), function(j) nay[id == j, , drop = FALSE])
   if (object$center){
      Xlist <- mapply(function(X, v){
         scale(X, X[v==0, ], scale = FALSE)
      }, Xlist, vlist, SIMPLIFY = FALSE)
      X <- do.call("rbind", Xlist)
   }
   Psi <- if (fixed == TRUE){
      diag(0, dim$q)
   } else {
      object$Psi
   }
   Xlist <- mapply(function(X, v) X[v!=0, , drop = FALSE], 
                   Xlist, vlist, SIMPLIFY = FALSE)
   ylist <- mapply(function(y, v) y[v!=0, , drop = FALSE], 
                   ylist, vlist, SIMPLIFY = FALSE)
   nalist <- mapply(function(na, v) na[v!=0], nalist, vlist, SIMPLIFY = FALSE)   
   gls <- glsfit(Xlist, Zlist = lapply(Xlist, function(z) z[, 1:dim$q, drop = FALSE]), 
                 ylist, Slist, nalist, Psi, onlycoef = FALSE)
   tlm <- summary(lm(gls$invtUy ~ 0 + gls$invtUX))
   tdata <- data.frame(gls$invtUy, gls$invtUX)
   names(tdata)[1] <- names(mf)[1]
   tdata$residuals <- tlm$residuals
   colnames(tdata) <- paste0("t", colnames(tdata))
   tdata$id <- id[v!=0]
   out <- list(tdata = tdata, 
               R2 = tlm$r.squared, R2adj = tlm$adj.r.squared,
               deviance = list(D = sum(tlm$residuals^2),
                               df = tlm$df[2],
                               p = pchisq(sum(tlm$residuals^2), tlm$df[2], lower.tail = F))
   )
   class(out) <- "gof.dosresmeta"
   return(out)
}


#' @rdname gof
#' @method print gof.dosresmeta
#' @export

print.gof.dosresmeta <- function(x, digits = 3, ...){
   signif <- symnum(x$deviance$p, corr = FALSE, na = FALSE, 
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    symbols = c("***", "**", "*", ".", " "))
   D <- formatC(x$deviance$D, digits = digits, format = "f")
   pvalue <- formatC(x$deviance$p, digits = digits, format = "f")
   R2 <- formatC(x$R2, digits = digits, format = "f")
   R2adj <- formatC(x$R2adj, digits = digits, format = "f")
   table <- cbind(D, df = x$deviance$df, `p-value` = pvalue, signif)
   colnames(table)[4] <- ""
   cat("Goodness-of-fit statistics:", "\n", sep = "")
   cat("\nDeviance test:", "\n")
   cat("D = ", D, " (df = ", x$deviance$df, "), p-value = ", pvalue, "\n", sep = "")
   cat("\n")
   cat("Coefficient of determination R-squared:", R2, "\n", sep = " ")
   cat("Adjusted R-squared:", R2adj, sep = " ")
   cat("\n")
}


#' Variance Partition Components for dosresmeta Objects
#' 
#' @name vpc
#' 
#' @description Computes the Variance Partition Components for dose-response meta-analysis.
#' 
#' @param object an object of class \code{dosresmeta} produced by \code{\link{dosresmeta}}.
#' 
#' @return A vector containing the variance partition components for each non-referent observation.
#' 
#' @examples
#' 
#' ## loading data
#' data("sim_os")
#' 
#' ## Quadratic (one-stage) dose-response model
#' quadr <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, type = type,
#'                     cases = cases, n = n, data = sim_os, proc = "1stage")
#'                     
#' ## Plot of the estimated vpc
#' plot(sim_os$dose[sim_os$se!=0], vpc(quadr), xlab = "dose")
#' lines(lowess(sim_os$dose[sim_os$se!=0], vpc(quadr)))
#'  
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references 
#' 
#' Goldstein H, Browne W, Rasbash J. Partitioning variation in multilevel models. 
#' Understanding Statistics: Statistical Issues in Psychology, Education, and the Social Sciences. 
#' 2002 Dec 2;1(4):223-31.
#' 
#' @export vpc

vpc <- function(object){
   v <- object$v
   id <- object$id
   Psi <- object$Psi
   Slist <- object$Slist
   vlist <- lapply(unique(id), function(j) cbind(v[id == j]))
   Z <- model.matrix(object$formula, data = object$model)[, 2:(object$dim$q+1), drop = FALSE]
   Zlist <- lapply(unique(id), function(j)
      Z[id == j, , drop = FALSE])
   if (object$center){
      Zlist <- mapply(function(Z, v){
         scale(Z, Z[v==0, ], scale = FALSE)
      }, Zlist, vlist, SIMPLIFY = FALSE)
      Z <- do.call("rbind", Zlist)
   }
   Zlist <- lapply(Zlist, function(z) z[-1, , drop = FALSE])
   vpclist <- mapply(function(Z, S){
      diag(Z %*% Psi %*% t(Z))/diag(S + Z %*% Psi %*% t(Z))
   }, Zlist, Slist, SIMPLIFY = FALSE)
   vpc <- do.call("c", vpclist)
   vpc
}


#' Grid with combinations of p for two-order fractional polynomials
#' 
#' @name fpgrid
#' 
#' @description Computes the different combinations of p usefull for evaluating two-order fractional polynomials.
#' 
#' @param p a numeric vector with the coefficient to be combined.
#' 
#' @return A data.frame with the different combinations of p.
#' 
#' @examples
#' 
#' grd <- fpgrid()
#' head(grd)
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references 
#' 
#' Royston, Patrick, and Douglas G. Altman. "Regression using fractional polynomials of 
#' continuous covariates: parsimonious parametric modelling." Applied Statistics (1994): 429-467.
#' 
#' @export fpgrid

fpgrid <- function(p = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)){
   if (length(p) < 2L){
      stop("p must contain at least 2 values")
   }
   p1 <- p
   p2 <- p   
   grid <- subset(expand.grid(p1 = p1, p2 = p2), p1 <= p2)
   grid <- grid[order(grid[, 1]), ]
   rownames(grid) <- seq(nrow(grid))
   grid
}

