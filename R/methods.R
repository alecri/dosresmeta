#' dosresmeta Object
#' 
#' @description An object returned by \code{dosresmeta} function, inheriting
#'   from class "\code{dosresmeta}", and representing a fitted dose-response
#'   (meta-analytical) model.
#'   
#' @name dosresmetaObject
#' @docType class
#'   
#' @return Objects of class "\code{dosresmeta}" are lists with defined
#' components. Dimensions of such components differs according to the choosen
#' procedure. For the one-stage analysis the dimensions refer to a one
#' dimensional outcome, \eqn{p} predictors and \eqn{m} studies used for fitting
#' the model. For the two-stage analysis the dimensions refer to \eqn{p}
#' outcome parameters, no predictor (only the intercept) and \eqn{m} studies.
#' The following components needs to be included in a legitimate mvmeta object: 
#' \tabular{ll}{ 
#' \code{coefficients} \tab a \eqn{p}-dimensional vector of the fixed-effects coefficients. \cr 
#' \code{vcov} \tab estimated \eqn{p \times p}{p x p} (co)variance matrix of the fixed-effects coefficients. \cr
#' \code{Psi} \tab for random-effects models, the estimated \eqn{p \times p}{p x p} between-study (co)variance matrix. \cr
#' \code{residuals} \tab a vector of residuals, that is observed minus fitted values. \cr
#' \code{fitted.values} \tab a vector of of fitted mean values. \cr
#' \code{df.residual} \tab the residual degrees of freedom. \cr
#' \code{rank} \tab the numeric rank of the fitted model. \cr
#' \code{logLik} \tab the (restricted) log-likelihood of the fitted model. \cr
#' \code{converged,niter} \tab for models with iterative estimation methods, logical scalar indicating if the algorithm eventually 
#' converged, and number of iterations. \cr
#' \code{control} \tab a list with the values of the control arguments used, as returned by \code{\link{dosresmeta.control}}. \cr
#' \code{method} \tab the estimation method. \cr
#' \code{dim} \tab list with the following scalar components: \eqn{m} (number of studies included in estimation, \eqn{k} (number of outcome parameters), 
#' \eqn{p} (number of predictors). \cr
#' \code{df} \tab list with the following scalar components: nall (number of observations used for estimation, excluding missing values), 
#' nobs (equal to nall, minus the number of fixed-effects coefficients in REML models), fixed (number of estimated fixed-effects coefficients), 
#' random (number of estimated (co)variance terms). \cr
#' \code{lab} \tab list with the following label vectors: \eqn{p} for the p predictors (including intercept). \cr
#' \code{model} \tab the model frame used for fitting. \cr
#' \code{call} \tab the function call. \cr
#' \code{formula} \tab the model supplied. \cr
#' \code{terms} \tab the \code{\link{terms}} object representing the fitted model. \cr
#' \code{proc} \tab the estimation procedure. \cr
#' \code{center} \tab if the desing matrix had been centered. \cr
#' \code{covariance} \tab how the (co)variance had been appproximated. \cr
#' \code{Slist} \tab list of approximated (co)variance matrices. \cr
#' \code{id} \tab identification vector of the studies. \cr
#' \code{v} \tab variances of the outcome values
#' }
#' 
#' @section Methods:
#' 
#' A number of methods functions are available for \code{\link{dosresmeta}} objects, most of them common to other regression functions. 
#' Specifically-written method functions are defined for predict (standard predictions). The qtest method performs the Cochran Q test for heterogeneity only for a two-stage analysis. 
#' Other methods have been produced for summary, logLik, coef, and vcov. Printing functions for the objects of classes defined above are also provided.
#' All the methods above are visible (exported from the namespace) and documented. In additions, several default method functions for regression are also 
#' applicable to objects of class "mvmeta", such as fitted, residuals, AIC, BIC and update, among others.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{dosresmeta-package}}, \code{\link{mvmetaObject}}

NULL


#' Extract Coefficients and (Co)Variance Matrix from dosresmeta Objects
#'
#' @description These method functions return the estimated fixed-effects coefficients and their (co)variance matrix for fitted 
#' dose-response meta-analytical models represented in objects of class "\code{dosresmeta}".
#'
#' @param object an object of class "\code{dosresmeta}".
#' @param format format of the returned object.
#' @param \dots further arguments passed to or from other methods.
#' @return For \code{coef}, a vector (default) or matrix with the estimated (fixed-effects) coefficients. 
#' For \code{vcov}, the (co)variance matrix of the estimated (fixed-effects) coefficients.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{coef}}, \code{\link{vcov}}, \code{\link{logLik.dosresmeta}}
#' 
#' @name coef.dosresmeta
#' 
#' @examples
#' ## Load data and run the model
#' data("alcohol_cvd")
#' model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                    se = se, cases = cases, n = n, data = alcohol_cvd) 
#'
#' ## Fixed-effect coefficients
#' coef(model)
#'
#' ## Fixed-effect (co)variance matrix
#' vcov(model)
#'
#'
#' @method coef dosresmeta  
#' @export


coef.dosresmeta <- function (object, format = c("vector", "matrix"), ...){
   coef <- object$coefficients
   format <- match.arg(format, c("vector", "matrix"))
   if (format == "matrix" || is.vector(coef)) 
      return(coef)
   names <- paste(rep(colnames(coef), each = nrow(coef)), 
                  rep(rownames(coef), ncol(coef)), sep = "")
   coef <- as.numeric(coef)
   names(coef) <- names
   return(coef)
}

#' @rdname coef.dosresmeta
#' @method vcov dosresmeta
#' @export

vcov.dosresmeta <- function(object, ...){
   return(object$vcov)
}


#' Extract Log-Likelihood from dosresmeta Objects
#'
#' This method function returns the log-likelihood for fitted dose-response models represented in objects of class "\code{dosresmeta}".
#' 
#' @param object an object of class "\code{dosresmeta}"
#' @param \dots further arguments passed to or from other methods.
#' @return A numeric scalar of class \code{"logLik"}.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' @seealso \code{\link{dosresmeta}}, \code{\link{dosresmeta-package}}, \code{\link{logLik}}
#' 
#' @examples
#' data("alcohol_crc")
#' 
#' ## Dose-response model assuming linearity
#' lin <- dosresmeta(formula = logrr ~ dose, type = type, id = id, se = se, 
#'            cases = cases, n = peryears, data = alcohol_crc, proc = "1stage")
#' 
#' ## Log-likelihood
#' ll <- logLik(lin)
#' ll
#' attributes(ll)
#' 
#' ## AIC and BIC
#' AIC(ll)
#' BIC(ll)
#' 
#' @method logLik dosresmeta  
#' @export

logLik.dosresmeta <- function (object, ...){
   val <- object$logLik
   attributes(val) <- object$df
   class(val) <- "logLik"
   val
}


#' Summarizing dosresmeta Models
#' 
#' @description Print and summary method functions for dose-response models represented in objects of class "\code{dosresmeta}".
#' 
#' @param object an object of class \code{dosresmeta} produced by \code{\link{dosresmeta}}.
#' @param x an object of class \code{dosresmeta} or \code{summary.dosresmeta} produced by \code{\link{dosresmeta}} or \code{summary.dosresmeta}, respectively.
#' @param ci.level the confidence level used for defining the confidence intervals for the estimates of the (fixed-effects) coefficients.
#' @param digits an integer specifying the number of digits to which printed results must be rounded.
#' @param \dots further arguments passed to or from other methods.
#'
#' @details the \code{print} method for class \code{dosresmeta} only returns basic information of the fitted model, namely the call, 
#' estimated (fixed-effects) coefficients, dimensions, and fit statistics (log-likelihood, AIC, BIC).
#' 
#' The \code{summary} method function computes additional statistics and tests, and produces a list object of class \code{summary.dosresmeta}. 
#' The \code{print} method function for this class, depending on the number of studies included in the analysis, shows additional information, 
#' such as tables reporting the estimates for the fixed and random-effects parts of the model, Chi-square test for model significance, 
#' Cochran Q test for heterogeneity and I-square. 
#' 
#' @return The \code{summary} method function for \code{dosresmeta} objects produces a list of class "\code{summary.dosresmeta}".
#' The components of the lists are some of those stored in the related \code{dosresmeta} object, plus the following:
#' \tabular{ll}{ 
#' \code{AIC} \tab the value of the Akaike information criterion for the fitted \code{dosresmeta} model, obtained through a call to \code{\link{AIC}}. \cr 
#' \code{BIC} \tab the value of the Bayesian information criterion for the fitted \code{dosresmeta} model, obtained through a call to \code{\link{BIC}} \cr
#' \code{corFixed} \tab the \eqn{p \times p}{p x p} correlation matrix of the fixed-effects coefficients, 
#' obtained from the (co)variance matrix \code{\link{vcov}} \cr
#' \code{corRandom} \tab the \eqn{p \times p}{p x p} correlation matrix of the random effects, obtained from the between-study (co)variance matrix \eqn{\Psi} \cr
#' \code{qstat} \tab results from the Cochran Q test for heterogeneity. \cr
#' \code{ci.level} \tab the confidence level used for defining the confidence intervals for the estimates of the fixed-effects coefficients. \cr
#' \code{chisq} \tab overall test similar to anova. \cr
#' }
#' As usual, the \code{print} method functions for classes "\code{dosresmeta}" and "\code{summary.dosresmeta}" do not return any value.
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{summary}}
#' 
#' @examples
#' ## Load data and run the model
#' data("alcohol_cvd")
#' model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                     se = se, cases = cases, n = n, data = alcohol_cvd) 
#' ## Defult print
#' model
#' ## Specify digits
#' print(model, digit = 2)
#' ## summary with 90th confidence intervals
#' summary(model, ci.level = .8)
#'
#' @rdname summary.dosresmeta
#' @method print dosresmeta
#' @export

print.dosresmeta <- function (x, digits = 4, ...){
   cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
       "\n\n", sep = "")
   cat("Fixed-effects coefficients:", "\n", sep = "")
   table <- formatC(x$coefficients, digits = digits, format = "f")
   print(table, quote = FALSE, right = TRUE, print.gap = 2)
   cat("\n")
   cat(x$dim$m, ifelse(x$dim$m == 1L, " study ", " studies, "), x$df$nall, 
       " values, ", x$df$fixed, " fixed and ", x$df$random, 
       " random-effects parameters", "\n", sep = "")
   table <- c(x$logLik, AIC(x), BIC(x))
   names(table) <- c("logLik", "AIC", "BIC")
   table <- formatC(table, digits = digits, format = "f")
   print(table, quote = FALSE, right = TRUE, print.gap = 2)
   cat("\n")
}

#' @rdname summary.dosresmeta
#' @method summary dosresmeta
#' @export

summary.dosresmeta <- function(object, ci.level = 0.95, ...){
   if (ci.level <= 0 || ci.level >= 1) 
      stop("'ci.level' must be within 0 and 1")
   coef <- object$coefficients
   vcov <- object$vcov
   dim <- object$dim
   Psi <- object$Psi
   lab <- object$lab
   coef <- as.numeric(coef)
   coef.se <- sqrt(diag(vcov))
   zval <- coef/coef.se
   zvalci <- qnorm((1 - ci.level)/2, lower.tail = FALSE)
   pvalue <- 2 * (1 - pnorm(abs(zval)))
   ci.lb <- coef - zvalci * coef.se
   ci.ub <- coef + zvalci * coef.se
   cilab <- paste(signif(ci.level, 2) * 100, "%ci.", c("lb", 
                                                       "ub"), sep = "")
   Wald.test <- try(waldtest(Sigma = vcov(object), b = coef(object),
                           Terms = 1:nrow(vcov))$chitest)
   tabfixed <- cbind(coef, coef.se, zval, pvalue, ci.lb, ci.ub)
   #dimnames(tabfixed) <- list(colnames(object$lab), 
   #                           c("Estimate", "Std. Error", "z", "Pr(>|z|)", cilab))
   dimnames(tabfixed)[2] <- list(c("Estimate", "Std. Error", "z", "Pr(>|z|)", cilab))
   corFixed <- vcov/outer(coef.se, coef.se)
   corRandom <- if (object$method != "fixed"){
      ran.sd <- sqrt(diag(Psi))
      Psi/outer(ran.sd, ran.sd)
   } else {
      NULL
   }
   qstat <- unclass(qtest.dosresmeta(object))
   keep <- match(c("vcov", "Psi", "df.res", "rank", 
                   "logLik", "converged", "niter", "negeigen", "method", 
                   "dim", "df", "lab", "na.action", "call", "terms",
                   "covariance", "proc"), names(object), 
                 0L)
   out <- c(list(coefficients = tabfixed), object[keep], 
            list(AIC = AIC(object), BIC = BIC(object), corFixed = corFixed, 
                 corRandom = corRandom, Wald.test = Wald.test, qstat = qstat, 
                 ci.level = ci.level))
   class(out) <- "summary.dosresmeta"
   return(out)
}

#' @rdname summary.dosresmeta
#' @method print summary.dosresmeta
#' @export

print.summary.dosresmeta <- function(x, digits = 4, ...){
   methodname <- c("reml", "ml", "fixed")
   methodlabel <- c("REML", "ML", "Fixed")
   covariancename <- c("h", "gl", "md", "smd", "user", "indep")
   covariancelabel <- c("Hamling", "Greenland & Longnecker", 
                        "Mean Differences", "Standardized Mean Differences",
                        "User defined", "Independent")
   cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
       "\n\n", sep = "")
   cat(ifelse(x$proc == "1stage", "One-stage ", 
              "Two-stage "),
       ifelse(x$method == "fixed", "fixed", 
                               "random"), "-effects meta-", 
       ifelse(x$dim$p > x$dim$q, "regression", "analysis"), 
       "\n", sep = "")
   if (x$method != "fixed"){
      cat("Estimation method: ", methodlabel[
         which(x$method == methodname)], "\n", sep = "")
   }
   if (x$covariance != "user"){
      cat("Covariance approximation: ", covariancelabel[
         which(x$covariance == covariancename)], "\n", sep = "")
   }
   chi2 <- formatC(x$Wald.test["chi2"], digits = digits, format = "f")
   pchi2 <- formatC(x$Wald.test["P"], digits = digits, format = "f")
   cat("\n")
   cat("Chi2 model: X2 = ", chi2, " (df = ", x$Wald.test["df"], "), p-value = ", 
       pchi2, "\n", sep = "")
   cat("\n")
   cat("Fixed-effects coefficients", "\n", sep = "")
   signif <- symnum(x$coefficients[, "Pr(>|z|)"], corr = FALSE, 
                    na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    symbols = c("***", "**", "*", ".", " "))
   tabletot <- formatC(x$coefficients, digits = digits, format = "f")
   tabletot <- cbind(tabletot, signif)
   colnames(tabletot)[7] <- ""
   #rownames(tabletot) <- x$lab$p
   print(tabletot, quote = FALSE, right = TRUE, print.gap = 2)
   cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
   k <- x$dim$q #ifelse(x$proc == "2stages", x$dim$k, x$dim$p)
   if (!x$method == "fixed"){
      cat("Between-study random-effects (co)variance components", 
          "\n", sep = "")
      sd <- formatC(cbind(`Std. Dev` = sqrt(diag(x$Psi))), 
                    digits = digits, format = "f")
      if (k == 1L) rownames(sd) <- ""
      if (k > 1L) {
         corRan <- x$corRandom
         corRan[upper.tri(x$corRan, diag = TRUE)] <- NA
         dimnames(corRan) <- NULL
         corRan <- format(corRan[-1, -ncol(corRan), drop = FALSE], 
                          digits = digits, format = "f")
         corRan <- rbind(x$lab$q[-k], corRan)
         colnames(corRan) <- c("Corr", rep("", k - 2))
         corRan[grep("NA", corRan)] <- ""
      } else {
         corRan <- NULL
      }
      print(cbind(sd, corRan), quote = FALSE, right = TRUE, 
            na.print = "", print.gap = 2)
      if (!is.null(x$negeigen) && x$negeigen > 0) {
         cat("(Note: Truncated estimate - ", x$negeigen, 
             " negative eigenvalues set to 0)", 
             "\n", sep = "")
      }
      cat("\n")
   }
   if (x$proc != "1stage"){
      Q <- formatC(x$qstat$Q, digits = digits, format = "f")
      pvalue <- formatC(x$qstat$pvalue, digits = digits, format = "f")
      i2 <- formatC(pmax((x$qstat$Q - x$qstat$df)/x$qstat$Q * 100, 
                         0), digits = 1, format = "f")
      cat("Univariate ", "Cochran Q-test for ", "residual ", 
          "heterogeneity:", "\n", sep = "")
      cat("Q = ", Q[1], " (df = ", x$qstat$df[1], "), p-value = ", 
          pvalue[1], "\n", sep = "")
      cat("I-square statistic = ", i2[1], "%", "\n\n", sep = "")
   }
   cat(x$dim$m, ifelse(x$dim$m == 1L, " study ", " studies, "), x$df$nall, 
       " values, ", x$df$fixed, " fixed and ", x$df$random, 
       " random-effects parameters", "\n", sep = "")
   if (x$method != "mm"){
      table <- c(x$logLik, x$AIC, x$BIC)
      names(table) <- c("logLik", "AIC", "BIC")
      table <- formatC(table, digits = digits, format = "f")
      print(table, quote = FALSE, right = TRUE, print.gap = 2)
   }
   cat("\n")
}


#' Cochran Q Test of Heterogeneity for dosresmeta Models
#' 
#' @description This method function performs a Cochran Q test of (residual) heterogeneity on fitted dose-response meta-analytical models 
#' represented in objects of class "\code{doseremeta}". It is implemented only for a two-stage approach and will return \code{NULL} otherwise.
#' 
#' @param object objects of classe "\code{dosresmeta}".
#' @param x an object of class "\code{qtest.dosresmeta}".
#' @param digits an integer specifying the number of digits to which printed results must be rounded.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @details In (multivariate) dose-response meta-analytical models, the test assesses the null hypothesis that the variability in the 
#' (multivariate) distribution of the outcomes is explained only in terms of estimation error in each study, measured by the within-study (co)variance matrices 
#' stored in the component \code{Slist} of \code{dosresmeta} objects. This is equal to test the hypothesis that the between-study (co)variance matrix is 
#' a zero matrix, and there is no random deviation in study-specific estimates.
#' 
#' @method qtest dosresmeta
#' @export

qtest.dosresmeta <- function (object, ...){
   if (object$proc == "1stage") return(NULL)
   mf <- model.frame(object)
   id <- object$id
   bi <- object$bi
   Si <- object$Si
   v <- object$v
   Xstar <- model.matrix(object$mod, do.call("rbind", 
                                             lapply(lapply(unique(id), function(j) 
                                                mf[v!=0 &id == j, , drop = FALSE]), head, 1)))
   nabi <- is.na(bi)
   bilist <- lapply(seq(object$dim$m), function(i) bi[i, ][!nabi[i, ]])
   Xstarlist <- lapply(seq(object$dim$m), function(i) 
      diag(1, object$dim$q)[!nabi[i, ], , drop = FALSE] %x% Xstar[i, , drop = FALSE])
   nabilist <- lapply(seq(object$dim$m), function(i) nabi[i, ])
   Psi <- diag(0, object$dim$q)
   Zstarlist <- lapply(seq(object$dim$m), function(j)
      diag(object$dim$q))
   gls <- glsfit(Xstarlist, Zstarlist, bilist, Si, nabilist, Psi, onlycoef = FALSE)
   Q <- drop(crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
   df <- with(object$df, nall - fixed)
   if (object$dim$q > 1L){
      Q <- c(Q, colSums(do.call("rbind", mapply(function(y, S, X, na) {
         comp <- rep(0, object$dim$q)
         comp[!na] <- as.vector((y - X %*% gls$coef)^2/diag(S))
         return(comp)
      }, bilist, Si, Xstarlist, nabilist, SIMPLIFY = FALSE))))
      df <- c(df, colSums(!nabi, na.rm = TRUE) - object$dim$q)
   }
   pvalue <- sapply(seq(length(Q)), function(i) 1 - pchisq(Q[i], df[i]))
   names(Q) <- names(df) <- names(pvalue) <- c(".all", if(object$dim$q > 1L) object$lab$q)
   qstat <- list(Q = Q, df = df, pvalue = pvalue, residual = object$dim$q)
   class(qstat) <- "qtest.dosresmeta"
   qstat
}

# qtest.dosresmeta <- function(object, ...){ 
#    mf <- model.frame(object)
#    mfmX <- object$model
#    y <- as.matrix(model.response(mfmX, "numeric"))
#    id <- object$id
#    nay <- is.na(y)
#    X <- model.matrix(attr(mfmX, "terms"), data = mfmX)
#    X <- X[, -grep("(Intercept)", colnames(X)), drop = FALSE]
#    Slist <- object$Slist
#    v <- object$v
#    dim <- object$dim
#    vlist <- lapply(unique(id), function(j) cbind(v[id == j]))
#    ylist <- lapply(unique(id), function(j) 
#       y[id == j, , drop = FALSE][!nay[j, ], , drop = FALSE])
#    Xlist <- lapply(unique(id),
#                    function(j) 
#                       X[id == j, , drop = FALSE][!nay[j, ], , drop = FALSE])
#    nalist <- split(nay, id)
#    if (object$center){
#       Xlist <- mapply(function(X, v){
#          scale(X, X[v==0, ], scale = FALSE)
#       }, Xlist, vlist, SIMPLIFY = FALSE)
#       X <- do.call("rbind", Xlist)
#    }
#    Psi <- diag(0, dim$q)
#    Xlist <- mapply(function(X, v) X[v!=0, , drop = FALSE], 
#                    Xlist, vlist, SIMPLIFY = FALSE)
#    ylist <- mapply(function(y, v) y[v!=0, , drop = FALSE], 
#                    ylist, vlist, SIMPLIFY = FALSE)
#    nalist <- mapply(function(na, v) na[v!=0], nalist, vlist, SIMPLIFY = FALSE)   
#    if (object$proc == "2stage"){
#       twoStage <- mapply(function(X, y, S, na)
#          dosresmeta.fixed(X, X, y, S, na, q = dim$q, 
#                           nall = sum(!nay), object$control),
#          lapply(Xlist, function(x) list(x[, 1:dim$q, drop = FALSE])),
#          lapply(ylist, function(x) list(x)), 
#          lapply(Slist, function(x) list(x)),
#          lapply(nalist, function(x) list(x)),
#          SIMPLIFY = FALSE)
#       y <- do.call("rbind", lapply(twoStage, function(x) x$coefficients))
#       dimnames(y) <- list(unique(id), object$lab$q)
#       ylist <- lapply(seq(dim$m), function(i) y[i, ][!nay[i, ]])
#       Slist <- lapply(twoStage, function(x) x$vcov)
#       Xlist <- lapply(seq(dim$m), 
#                       function(i) diag(1, dim$q) %x% 
#                          model.matrix(object$mod, 
#                                       data = do.call("rbind", lapply(split(mfmX, id), head, 1)))[i, , drop = FALSE]
#       )
#       Psi <- diag(0, dim$q)
#       nay <- is.na(y)
#       nalist <- split(nay, unique(id))      
#    }
#    gls <- glsfit(Xlist, Zlist = lapply(Xlist, function(z) z[, 1:dim$q, drop = FALSE]), 
#                  ylist, Slist, nalist, Psi, onlycoef = FALSE)
#    Q <- drop(crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
#    df <- with(object$df, nall - fixed)
#    pvalue <- 1 - pchisq(Q, df)
#    qstat <- list(Q = Q, df = df, pvalue = pvalue)
#    class(qstat) <- "qtest.dosresmeta"
#    qstat
# }

#' @rdname qtest.dosresmeta
#' @method print qtest.dosresmeta
#' @export

print.qtest.dosresmeta <- function (x, digits = 3, ...){
   signif <- symnum(x$pvalue, corr = FALSE, na = FALSE, 
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    symbols = c("***", "**", "*", ".", " "))
   Q <- formatC(x$Q, digits = digits, format = "f")
   pvalue <- formatC(x$pvalue, digits = digits, format = "f")
   table <- cbind(Q, df = x$df, `p-value` = pvalue, signif)
   colnames(table)[4] <- ""
   cat("Cochran Q-test for ", "heterogeneity", "\n", sep = "")
   cat("\nOverall test:", "\n")
   cat("Q = ", Q[1], " (df = ", x$df[1], "), p-value = ", pvalue[1], 
       "\n", sep = "")
   cat("\n")
}


#' Predicted Values from dosresmeta Models 
#'
#' @description This method function computes predictions from fitted dose-response models 
#' represented in objects of class "\code{dosresmeta}", optionally for a new set of exposure levels. 
#' Predictions are optionally accompanied by confidence intervals and/or standard errors for the predictions.
#' 
#' @param object an object of class \code{dosreseta}.
#' @param newdata an optional data frame or matrix in which to look for variables values with which to predict from dose-response models.
#' @param xref an optional scalar to indicate which levels should serve as referent for the predicted relative risks. See details.
#' @param expo logical switch indicating if the prediction should be on the exponential scale.
#' @param xref_vec an optional numeric to indicate the referent (vector) for the predicted relative risks. See details.
#' @param se.incl logical switch indicating if standard errors need to be included.
#' @param ci.incl logical switch indicating if confidence intervals need to be included.
#' @param ci.level a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.
#' @param xref_pos an optional scalar to indicate the position of the referent for the predicted relative risks. See details.
#' @param order logical to indicate if the predictions need to be sorted by exposure levels.
#' @param delta an optional scalar to specify to predict the linear trend related to that increase.
#' @param \dots further arguments passed to or from other methods.
#' @details The method function \code{predict} produces predicted values from \code{dosresmeta} objects. When more than one study is included in the analysis,
#' estimated predictions are only based on the fixed part of the model.
#' 
#' If \code{newdata} is omitted, the predictions are based on the data used for the fit. If \code{xref} is provided, it must be equal to one of the modeled values. 
#' If not
#' provided, the minimum modeled referent value will be used as referent for the predicted relative risks
#' 
#' If \code{newdata} is specified, it should include all the variables used to model the dose-response relation. Again, if specified, \code{xref} must be equal to one
#' of the value in the newdata. If omitted, the minimum value for the newdara will be used as referent.
#' 
#' Only for the linear trend it is possible to specify the predicted increase of risk correspongind to an increase equal to \code{delta} argument.
#' 
#' By default (\code{order = TRUE}), the predictions are sorted by exposure levels to facilitate understanding and possible graphical
#' presentation of the results.
#' 
#' @return The results are returned structured in a data frame.
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' @seealso \code{\link{dosresmeta}}, \code{\link{predict}}
#' 
#' @examples
#' ## Load data and run the linear and quadratic models
#' data("alcohol_cvd")
#' lin <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                   se = se, cases = cases, n = n, data = alcohol_cvd) 
#' quadr <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                     se = se, cases = cases, n = n, data = alcohol_cvd) 
#' 
#' ## Predicted linear trend (on RR scale)
#' predict(lin, delta = 12, expo = TRUE)
#' 
#' ## Predicted modeled data from quadratic model (on RR scale)
#' predict(quadr, expo = TRUE)
#' 
#' ## Plot predicted dose-response relation
#' with(predict(quadr, order = TRUE, expo = TRUE), {
#'    plot(dose, pred, log = "y", type = "l",
#'         xlim = c(0, 45), ylim = c(.4, 2))
#'    lines(dose,  ci.lb, lty = 2)
#'    lines(dose, ci.ub, lty = 2)
#'    rug(dose, quiet = TRUE)
#' })
#' 
#' ## Prediction for new values from quadratic model (on RR scale)
#' newdata <- data.frame(dose = seq(0, 50, 5))
#' predict(quadr, newdata, expo = TRUE)
#' 
#' ## Loading data
#' data("ari")
#' 
## Quadratic model for standardized mean differences
#' mod <- dosresmeta(y ~ dose + I(dose^2), id = id, sd = sd,
#'  n = n, data = ari, covariance = "smd")
#' 
#' ## Smoothed plot
#' newdata <- data.frame(dose = seq(0, 30, 1))
#' with(predict(mod, newdata), {
#'    plot(dose, pred, type = "l", ylim = c(0, .6))
#'    lines(dose,  ci.lb, lty = 2)
#'    lines(dose, ci.ub, lty = 2)
#'    rug(dose, quiet = TRUE)
#' })
#' 
#' @rdname predict.dosresmeta
#' @method predict dosresmeta
#' @export

predict.dosresmeta <- function(object, newdata, xref, expo = FALSE, xref_vec,
                               ci.incl = TRUE, se.incl = FALSE, xref_pos = 1,
                               delta, order = FALSE, ci.level = 0.95, ...){
   if (!missing(delta)){
      if (object$dim$q > 1L)
         stop("'delta' option available only for linear trend")
      mf <- model.frame(object)
      xref_delta <- if (missing(xref_vec)){
         if (missing(xref)){
            0
         } else {
            xref
         }
      } else {
         xref_vec
      }
      mf[1, 2] <- xref_delta
      mf[2, 2] <- xref_delta + delta
      X <- model.matrix(attr(mf, "terms"), data = mf)[1:2, , drop = FALSE]
      X <- X[, -grep("Intercept", colnames(X)), drop = FALSE]
   } else if (missing(newdata) || is.null(newdata)){
      mf <- model.frame(object)
      X <- model.matrix(attr(mf, "terms"), data = mf)[, -1, drop = FALSE]
   }  else {
      ttnr <- delete.response(terms(object))
      mf <- model.frame(ttnr, newdata)
      X <- model.matrix(ttnr, mf, contrasts.arg = object$contrasts)[, -1, drop = FALSE]
   }
   xref <- if (missing(xref_vec)){
      if (missing(xref)){
         X[xref_pos, ]
      } else {
         X[X[, 1] == xref, , drop = FALSE][1, ]
      }
   } else {
      xref_vec
   }
   if (!missing(delta)){
      xref <- xref_delta
   }
   fit <- X
   if (object$center){
      X <- scale(X, xref, scale = FALSE)
   }
   pred <- tcrossprod(X, rbind(c(t(object$coefficients))))
   fit <- if (expo == T) {
      cbind(fit, exp(pred))
   } else {
      cbind(fit, pred)
   }
   colnames(fit) <- c(colnames(X), "pred")
   zvalci <- qnorm((1 - ci.level)/2, lower.tail = FALSE)
   if (ci.incl == T){
      se <- sqrt(diag(X %*% tcrossprod(object$vcov, X)))
      ci.lb <- pred - zvalci * se
      ci.ub <- pred + zvalci * se
      fit <- if (expo == T){
         cbind(fit, exp(ci.lb), exp(ci.ub))
      } else {
         cbind(fit, ci.lb, ci.ub)
      }
      colnames(fit) <- c(colnames(X), "pred", "ci.lb", "ci.ub")
   }
   if (se.incl == T){
      fit <- cbind(fit, se = se)
   }
   fit <- as.data.frame(fit)
   if (order == T){
      fit <- fit[order(fit[, 1]), ]      
   }
   if (!missing(delta)){
      fit[1] <- delta
      fit <- fit[fit["ci.lb"] != fit["ci.ub"], ]
      colnames(fit)[1] <- "delta"
      row.names(fit) <- ""
   }
   fit
}

#' Best Linear Unbiased Predictions from dosresmeta Models
#' 
#' @description This method function computes (empirical) best linear unbiased predictions 
#' from fitted dose-response meta-analytical models represented in objects of class "dosresemta". 
#' 
#' @param object objects of classe "\code{dosresmeta}".
#' @param \dots further arguments passed to or from other methods.
#' 
#' @details The method function blup produces (empirical) best linear unbiased predictions from dosresmeta objects. 
#' Predictions are expressed in terms of study-specific deviations as random effects.
#' Predicted random effects from blup are a shrunk version of study-specific realizations, where study-specific predictions borrow 
#' strength from the assumption of an underlying distribution in a (usually hypothetical) population of studies.
#' Blup are not avaialable for fixed-effects models since the are meaningless in that context.
#' 
#' @examples
#' ## Load data and run the linear and quadratic models
#' data("alcohol_cvd")
#' lin <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                   se = se, cases = cases, n = n, data = alcohol_cvd) 
#' quadr <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                     se = se, cases = cases, n = n, data = alcohol_cvd) 
#' 
#' ## blup prediction for the previous models
#' blup(lin)
#' blup(quadr)
#' 
#' @method blup dosresmeta
#' @export 

blup.dosresmeta <- function(object, ...){
   if (object$method == "fixed") stop("blup predictions are available only for random-effects model")
   v <- object$v
   id <- object$id
   Psi <- object$Psi
   vlist <- lapply(unique(id), function(j) cbind(v[id == j]))
   res <- object$residuals
   if (object$proc == "2stage"){
      reslist <- lapply(seq(object$dim$m), function(j)
         t(cbind(res)[j, , drop = FALSE]))
      Zlist <- lapply(seq(object$dim$m), function(j)
         diag(object$dim$q))
      Vlist <- mapply(function(S, Z) 
         (S + Z %*% object$Psi %*% t(Z)), 
         object$Si, Zlist, SIMPLIFY = FALSE)
   } else if (object$proc == "1stage"){
      reslist <- lapply(unique(id), function(j) 
         res[id[v!=0] == j, , drop = FALSE])
      #Z <- object$model[, -1, drop = FALSE]
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
      Vlist <- mapply(function(S, Z) 
         (S + Z %*% object$Psi %*% t(Z)), 
         object$Slist, Zlist, SIMPLIFY = FALSE)
   }
   Ulist <- lapply(Vlist, chol)
   invUlist <- lapply(Ulist, function(U) backsolve(U, diag(ncol(U))))
   bluplist <- mapply(function(z, invU, res) {
      blup <- drop(Psi %*% t(z) %*% (tcrossprod(invU)) %*% res)
      names(blup) <- object$lab$q
      return(blup)
   }, Zlist, invUlist, reslist, SIMPLIFY = FALSE)
   do.call("rbind", bluplist)
}