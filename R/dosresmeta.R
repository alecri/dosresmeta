#' Multivariate Dose-Response Meta-Analysis
#' 
#' @description The function \code{dosresmeta} estimates a dose-response curve from either single or multiple summarized dose-response data, taking into account 
#' the correlation among observations and heterogeneity across studies. The function \code{dosresmeta.fit} is a wrapper for actual fitting functions based on 
#' different estimation methods, usually called internally. See \code{\link{dosresmeta-package}} for an overview.
#' 
#' @param formula an object of class "\code{\link{formula}}" offering a symbolic description of the dose-response functional relation. Terms in the formula can be
#' provided in the \code{data} below.
   #' @param mod an object of class "\code{\link{formula}}" offering a symbolic description of the meta-regression model (by default \code{mod = ~ 1}). Terms in the formula can be
#' provided in the \code{data} below.
#' @param X processed design matrix of fixed effects.
#' @param Z processed design matrix of random effects.
#' @param y processed outcome vector.
#' @param Slist list of approximated or given (co)variance matrices.
#' @param id an vector to specify the id variable for the studies included in the analysis. Optional if estimating a dose-response model from a single study.
#' @param v a vector to specify the variances of the reported outcome. Alternatively the user can provide the standard error in the \code{se} argument, 
#' or only for log relative risks, the confidence interval in the \code{lb} and \code{ub} arguments.
#' @param type an optional vector (or a string) required when the outcome is log relative risks. It specifies the study-specific design. 
#' The values for case-control, incidence-rate, and cumulative incidence data are \code{cc}, \code{ir}, and \code{ci}, respectively.
#' @param cases a vector to specify the number of cases for each exposure level. Required to reconstruct the (co)variance matrix for log relative risks.
#' @param n a vector to specify the total number of subjects for each exposure level. Required to reconstruct the (co)variance matrix for log relative risks.
#' For incidence-rate data \code{n} indicates the amount of person-time for each exposure level.
#' @param sd a vector to specify the standard deviation. Required to reconstruct the (co)variance matrix for differences and standardized mean differences.
#' @param data a data frame (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the previous arguments.
#' @param intercept a logical value to specify if an intercept term needs to be included in the model. See details.
#' @param center a logical value to specify if the design matrix need to be center at the referent ones. See details.
#' @param se an optional vector to specify the standard error of the reported log relative risks; needed if \code{v} is not provided.
#' @param lb an optional vector to specify the lower bound of the confidence interval for the reported relative risks; needed if \code{v} and \code{se} are not provided.
#' @param ub an optional vector to specify the upper bound of the confidence interval for the reported relative risks; needed if \code{v} and \code{se} are not provided.
#' @param covariance method to approximate the (co)variance matrix of the outcome. Options are "\code{gl}" for the method proposed by Greenland and Longnecker (default)
#' , "\code{h}" for the method proposed by Hamling, "\code{md}" for mean differences, "\code{smd}" for standardized mean differences, and "\code{user}" 
#' if provided by the user.
#' @param method.smd character string indicating the method to be used. Options are \code{cohens}, \code{hedges}, and \code{glass}. Required only if \code{covariance}
#' equal "\code{smd}".
#' @param proc "\code{2stage}" (default) or "\code{1stage}" procedure. See \code{\link{dosresmeta-package}} for an overview.
#' @param method method used to estimate the (pooled) dose-response relation: "\code{fixed}" for fixed-effects models, "\code{ml}" or "\code{reml}" 
#' for random-effects models fitted through (restricted) maximum likelihood, and "\code{mm}" for random-effects models fitted through method of moments (currently
#' available only for the two stages procedure). 
#' @param control list of parameters for controlling the fitting process. These are passed to \code{\link{dosresmeta.control}} by \code{\link{dosresmeta.fit}}
#' to replace otherwise selected default values.
#' 
#' @details The function defines all the elements required to estimate a dose-response association taking into account the correlation among the observations. 
#' If the (co)variance matrix is not provided then it is approximated depending of the type of outcome specified through the \code{covariance} argument.
#' The dose-response model is specified in the \code{formula}. Typically the outcome is expressed as a contrast to a reference exposure level, so that the model 
#' does not have an intercept and the values in the design matrix need to be centered at the referent values, as described by Qin Liu et al, 2009. 
#' This is internally done, respectively, when \code{intercept = FALSE} and \code{center = TRUE} (default values).
#' 
#' The function calls the wrapper \code{dosresmeta.fit} to perform the actual fitting. The latter prepares the data and calls specific fitting functions, 
#' depending on the chosen procedure and method. For the two stages procedure, the second part of the analysis is performed using the function \code{\link{mvmeta.fit}} 
#' from the \code{\link{mvmeta}} package. Different estimator are implemented in the package. The estimation options available are
#' \itemize{
#' \item Fixed-effects
#' \item Maximum likelihood (ML)
#' \item Restricted maximum likelihood (REML)
#' \item Method of moments (currently available only for the two stage procedure)
#' }
#' 
#' The fitting procedure can be controlled through the additional terms specified in control, which are passed to the function \code{\link{dosresmeta.control}}.
#' 
#' @return The \code{dosresmeta} function typically returns a list of object of class \code{dosresmeta} representing the meta-analytical model fit, 
#' as described in \code{\link{dosresmetaObject}}.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. 
#' American journal of epidemiology, 135(11), 1301-1309.
#' 
#' Orsini, N., Bellocco, R.,  Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.
#'
#' Liu, Q., Cook, N. R., Bergstrom, A., Hsieh, C. C. (2009). A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear 
#' dose-response data. Computational Statistics & Data Analysis, 53(12), 4157-4167. 
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. 
#' Statistics in Medicine, 31(29), 3821-3839.
#' 
#' @seealso \code{\link{dosresmeta-package}}, \code{\link{mvmeta}}, \code{\link{covar.logrr}}, \code{\link{covar.smd}}
#' 
#' @export dosresmeta
#' 
#' @examples
#' ## First example: Single case-control study
#' ## Linear trend estimation
#' data("cc_ex")
#' 
#' ## Fitting the model
#' mod1 <- dosresmeta(formula = logrr ~ dose, type = "cc", cases = case,
#'                    n = n, lb = lb, ub = ub, data= cc_ex)
#' summary(mod1)
#' ## Results
#' predict(mod1, delta = 1, expo = TRUE)
#' 
#' 
#' ## Second example: Multiple studies
#' ## Linear and quadratic trend using random-effects meta-analysis
#' data("alcohol_cvd")
#' 
#' ## Linear trend
#' lin <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                   se = se, cases = cases, n = n, data = alcohol_cvd)
#' summary(lin)
#' ## Predicted linear trend
#' predict(lin, delta = 1, expo = TRUE)
#' 
#' ## Non-linear (quadratic) trend
#' quadr <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                     se = se, cases = cases, n = n, data = alcohol_cvd)
#' summary(quadr)
#' 
#' ## Graphical results
#' with(predict(quadr, expo = TRUE, order = TRUE), {
#'    plot(dose, pred, log = "y", type = "l",
#'         xlim = c(0, 45), ylim = c(.4, 2))
#'    lines(dose,  ci.lb, lty = 2)
#'    lines(dose, ci.ub, lty = 2)
#'    rug(dose, quiet = TRUE)
#' })
#' 
#' 
#' ## Third example: Continous outcome (smd)
#' data("ari")
#' mod3 <- dosresmeta(formula = y ~ dose + I(dose^2), id = id,
#'                    sd = sd, n = n, covariance = "smd", data = ari)
#' summary(mod3)
#' 
#' ## Graphical results
#' newdata <- data.frame(dose = seq(0, 30, 1))
#' with(predict(mod3, newdata, order = TRUE), {
#'    plot(dose, pred, type = "l",
#'         ylim = c(0, .6))
#'    lines(dose,  ci.lb, lty = 2)
#'    lines(dose, ci.ub, lty = 2)
#'    rug(dose, quiet = TRUE)
#' })

dosresmeta <- function(formula, id, v, type,  cases, n, sd, data, 
                       mod = ~ 1, intercept = F, center = T, se, lb, ub, 
                       covariance = "gl", method = "reml", proc = "2stage", 
                       Slist, method.smd = "cohen", control = list()){
   if (missing(data)) 
      data <- NULL
   if (is.null(data)){
      data <- sys.frame(sys.parent())
   }
   else{
      if (!is.data.frame(data)){
         data <- data.frame(data)
      }
   }
   call <- match.call()
   mf <- match.call(expand.dots = FALSE)
   mf.id <- mf[[match("id", names(mf))]]
   ## ID doesn't need to be specified in case of single study analysis, so it is 
   ## set equal to one
   id <- if(is.null(mf.id)){
      as.factor(rep(1, nrow(data)))      
   } else {
      as.factor(eval(mf.id, data, enclos = sys.frame(sys.parent())))      
   }
   covariance <- match.arg(covariance, c("gl", "h", "md", "smd", "user", "indep"))
   proc <- match.arg(proc, c("1stage", "2stage"))
   method <- if (proc == "1stage"){
      match.arg(method, c("fixed", "ml", "reml"))
   } else {
      match.arg(method, c("fixed", "ml", "reml", "mm", "vc"))
   }
   ## there is no 2stage  (nor (re)ml) procedure for single study analysis
   if (length(unique(id)) == 1L){
      proc <- "1stage"
      method <- "fixed"
   }
   ## Evaluating arguments only if needed
   if (covariance %in% c("md", "smd")){
      if (missing(sd) | missing(n)){
         stop("Arguments md and smd are required when covariance equal to 'md' or 'smd'")
      }
      mf.sd <- mf[[match("sd", names(mf))]]
      mf.n <- mf[[match("n", names(mf))]]
      sd <- eval(mf.sd, data, enclos = sys.frame(sys.parent()))
      n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
      
   }
   if (covariance %in% c("gl", "h")){
      if (missing(cases) | missing(n) | missing(type)){
         stop("Arguments cases, n, and type are required when covariance equal to 'gl' or 'h'")
      }
      mf.cases <- mf[[match("cases", names(mf))]]
      mf.n <- mf[[match("n", names(mf))]]
      mf.type <- mf[[match("type", names(mf))]]      
      type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
      cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
      n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
      ## For single study analysis type can be specified by a string
      if (length(type) == 1L){
         type <- rep(as.vector(mf.type), length(id))         
      }
   }
   ## formulaMod specifies the model for the fixed-effects part
   formulaMod <- if (mod != ~ 1){
      ## problem if the user specify no intercept
      formula <- update(formula, . ~ . + 1)
      update(formula, paste(". ~ . +", paste(
         paste(formula[3], attr(terms(mod), "term.labels"), sep = ":"), 
         collapse = " + ")))
   } else {
      formula
   }
   ## it should be not possible otherwise (no intercept model)
   if (intercept == FALSE){
      formula <- update(formula, . ~ . + 0)   
   }
   mf.v <- mf[[match("v", names(mf))]]
   mf.se <- mf[[match("se", names(mf))]]
   mf.lb <- mf[[match("lb", names(mf))]]
   mf.ub <- mf[[match("ub", names(mf))]]
   if (covariance %in% c("gl", "h")){
      if (is.null(mf.v) & is.null(mf.se) & is.null(mf.lb) & is.null(mf.ub)){
         stop("One of the arguments v, se, or both lb and ub, needs to be specifed")
      }
   }
   ## probably to be simplified (less flexible)
   if (is.null(mf.v)){
      se <- eval(mf.se, data, enclos = sys.frame(sys.parent()))
      if (any(se < 0, na.rm = TRUE)){
         stop("Standard errors cannot be negative")
      }
      if (is.null(mf.se) && !(is.null(mf.lb) | is.null(mf.ub))){
         lb <- eval(mf.lb, data, enclos = sys.frame(sys.parent()))
         ub <- eval(mf.ub, data, enclos = sys.frame(sys.parent()))
         v <- ((log(ub) - log(lb))/(2 * qnorm(1 - .05/2)))^2
      }
      else {
         v <- se^2
      }
   }
   else {
      v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
   }
   v[is.na(v)] <- 0
   if (any(v < 0, na.rm = TRUE)){
      stop("Variances cannot be negative")
   }
   if (sum(v == 0) > 1L & length(unique(id)) == 1L){
      stop("It seems multiple studies are including without corresponding id argument. Please check.")
   }
   ## model frame for: fixed-effects (mfmX) and random-effects (mfmZ)
   mfmX <- model.frame(formulaMod, data)
   mfmZ <- model.frame(formula, data)
   terms <- attr(mfmX, "terms")
   y <- as.matrix(model.response(mfmX, "numeric"))
   nay <- is.na(y)
   ylist <- lapply(unique(id), function(j) y[id == j, ][!nay[j, ]])
   vlist <- lapply(unique(id), function(j) cbind(v[id == j]))
   nalist <- lapply(unique(id), function(j) nay[id == j])
   if (covariance %in% c("gl", "h")){
      caseslist <- lapply(unique(id), function(j) cases[id == j])
      nlist <- lapply(unique(id), function(j) n[id == j])
      typelist <- lapply(unique(id), function(j) type[id == j])
      Slist <- mapply(function(cases, n, y, v, type) 
         covar.logrr(cases, n, y, v, type, covariance = covariance), 
         caseslist, nlist, ylist, vlist, typelist, SIMPLIFY = FALSE)
   }
   if (covariance %in% c("md", "smd")){
      sdlist <- lapply(unique(id), function(j) sd[id == j])
      nlist <- lapply(unique(id), function(j) n[id == j])
      if (covariance == "smd") 
         method.smd <- match.arg(method.smd, c("cohens", "hedges", "glass"))
      covlist <- mapply(function(y, sd, n)
         covar.smd(y, sd, n, covariance, method.smd),
         ylist, sdlist, nlist, SIMPLIFY = FALSE)
      y <- do.call("rbind", lapply(covlist,  function(x) x$y))
      mfmX[, 1] <- y
      v <- cbind(unlist(lapply(covlist, function(x) x$v)))
      vlist <- lapply(unique(id), function(j) cbind(v[id == j]))
      Slist <- lapply(covlist, function(x) x$S)
   }
   if (covariance == "user"){
      if (!is.list(Slist)){
         Slist <- list(Slist)
      }
      if (any(sapply(Slist, function(s) nrow(s) != ncol(s)))){
         stop ("At least one of the covariance matrix is not squared")
      }
      if (any(do.call("c", Map(function(v, s){
         length(v) != (nrow(s) + 1)
      }, vlist, Slist)))){
         stop ("The dimension of Slist does not match the data")
      }
   }
   if (covariance == "indep"){
      Slist <- lapply(vlist, function(x) diag(x[x!=0], nrow = sum(x!=0)))
   }
   X <- model.matrix(attr(mfmX, "terms"), data = mfmX)[, -1, drop = FALSE]
   Xlist <- lapply(unique(id), function(j)
      X[id == j, , drop = FALSE][!nay[j, ], , drop = FALSE])
   Z <- model.matrix(attr(mfmZ, "terms"), data = mfmZ)
   Zlist <- lapply(unique(id), function(j)
      Z[id == j, , drop = FALSE][!nay[j, ], , drop = FALSE])
   if (center == T){
      Xlist <- mapply(function(X, v) scale(X, X[v==0, ], scale = FALSE), 
                      Xlist, vlist, SIMPLIFY = FALSE)
      X <- do.call("rbind", Xlist)
      Zlist <- mapply(function(Z, v) scale(Z, Z[v==0, ], scale = FALSE), 
                      Zlist, vlist, SIMPLIFY = FALSE)
      Z <- do.call("rbind", Zlist)
   }
   ## using only non-referent values
   fit <- dosresmeta.fit(X[v!=0, , drop = FALSE], Z[v!=0, , drop = FALSE], 
                         y[v!=0, , drop = FALSE], Slist, id[v!=0], method, 
                         control, proc, mod, v, data)
   fit$call <- call
   #fit$formula <- update(formulaMod, . ~ . + 1)
   fit$formula <- formulaMod
   fit$model <- mfmX
   fit$mod <- mod
   fit$terms <- terms
   fit$contrasts <- attr(X, "contrasts")
   fit$proc <- proc
   fit$center <- center
   fit$covariance <- covariance
   fit$Slist <- Slist
   fit$id <- id
   fit$v <- v
   class(fit) <- "dosresmeta"
   fit
}


#' @rdname dosresmeta
#' @export
dosresmeta.fit <- function (X, Z, y, Slist, id, method, control, 
                            proc, mod, v, data){
   control <- do.call("dosresmeta.control", control)
   m <- length(unique(id))
   k <- table(id)
   p <- ncol(X)
   q <- ncol(Z)
   nay <- is.na(y)
   n <- nrow(y)
   nall <- sum(!nay)
   ## names for coefficient
   np <- colnames(X)
   if (is.null(np))
      np <- paste("X", seq(p), sep = "")
   nq <- colnames(Z)
   if (is.null(nq))
      nq <- paste("Z", seq(q), sep = "")
   ny <- colnames(y)
   if (is.null(ny))
      ny <- "y"
   ylist <- lapply(unique(id), function(j) y[id == j, ][!nay[j, ]])
   Xlist <- lapply(unique(id), function(j) X[id == j, , drop = FALSE][
      !nay[j, ], , drop = FALSE])
   Zlist <- lapply(unique(id), function(j) Z[id == j, , drop = FALSE][
      !nay[j, ], , drop = FALSE])
   nalist <- lapply(unique(id), function(j) nay[id == j])
   if (proc == "2stage"){
      if (any(k < p)){
         stop("A two-stage approach requires that each study provides at least p non-referent obs (p is the number of columns of the design matrix X)")
      }
      ## obs: the argument of dosresmeta.fixed need to be list and since mapply
      ## requires additional list, I had to use lapply(. , . list(.))
      twoStage <- mapply(function(X, y, S, na)
         dosresmeta.fixed(X, X, y, S, na, q, nall, control),
         lapply(Zlist, function(z) list(z)),
         lapply(ylist, function(x) list(x)), 
         lapply(Slist, function(x) list(x)),
         lapply(nalist, function(x) list(x)),
         SIMPLIFY = FALSE)
      beta <- do.call("rbind", lapply(twoStage, function(x) x$coefficients))
      dimnames(beta) <- list(unique(id), nq)
      Sbeta <- do.call("list", lapply(twoStage, function(x) x$vcov))
      Xstar <- model.matrix(mod, do.call("rbind", 
                                         lapply(lapply(unique(id), function(j) data[v!=0, ][id == j, ]), head, 1)))
      fit <- mvmeta.fit(X = Xstar, y = beta, 
                        S = do.call("rbind", lapply(Sbeta, vechMat)), 
                        method = method)
      fit$bi <- beta
      fit$Si <- Sbeta
      fit[c("bscov", "xlevels", "S")] <- NULL
      names(fit$lab) <- c("q", "p")
      if (q == p) fit$lab$p <- nq
      fit
   } else {
      fun <- paste("dosresmeta", method, sep = ".")
      fit <- do.call(fun, list(Xlist = Xlist, Zlist = Zlist, ylist = ylist, 
                               Slist = Slist, nalist = nalist, q = q, 
                               nall = nall, control = control))
      if (!is.null(fit$converged) && !fit$converged){
         warning("convergence not reached after maximum number of iterations")
      }
      fit$method <- method
      fit$df <- list(nall = nall, nobs = nall - (method == "reml") * 
                        fit$rank, df = nall - fit$df.residual, fixed = fit$rank, 
                     random = ifelse(method == "fixed", 0, nall - fit$rank - 
                                        fit$df.residual))
      fit$lab <- list(p = np, q = nq)
      fit$fitted.values <- cbind(unlist(fit$fitted.values))
      if (method != "fixed") dimnames(fit$Psi) <- list(nq, nq)
      names(fit$coefficients) <- np
      rownames(fit$vcov) <- colnames(fit$vcov) <- np
      fit$residuals <- y - fit$fitted.values
      dimnames(fit$residuals)[2] <- dimnames(fit$fitted.values)[2] <- list(ny)
   }
   fit$dim <- list(m = m, p = p, q = q, k = k)
   fit
}


#' Ancillary Parameters for Controlling the Fit in dosresmeta Models
#' 
#' @description This internal function sets the parameter options used for fitting dose-response meta-analytical models, commonly to pre-specified default values. 
#' It is usually internally called by \code{\link{dosresmeta.fit}}.
#' 
#' @param optim list of parameters passed to the control argument of the function optim, which performs the quasi-Newton optimization in likelihood-based 
#' random-effects models. See \code{\link{optim}}.
#' @param showiter logical. If \code{TRUE}, the progress of iterative optimization is shown.
#' @param gr indicates if the gradient of the (re)ml likelihood should be provided. FALSE by default.
#' @param maxiter positive interger value. Maximum number of iterations in methods involving optimization procedures.
#' @param initPsi either a matrix or a vector of its lower triangular elements (with diagonal, taken by column) from which starting 
#' values of the parameters of the between-study (co)variance matrix are derived, used in the optimization procedure for likelihood-based random-effects models. 
#' If \code{NULL} (the default, and recommended), the starting value is created internally through an iterative generalized least square algorithm.
#' @param igls.iter number of iteration of the iterative generalized least square algorithm to be run in the hybrid optimization procedure 
#' of linkelihood-based models to provide the starting value.
#' @param reltol relative convergence tolerance in methods involving optimization procedures. The algorithm stops if it is unable to 
#' reduce the value by a factor of \eqn{reltol * (abs(val) + reltol)} at a step.
#' @param set.negeigen positive value. Value to which negative eigenvalues are to be set in estimators where such method is used
#' to force positive semi-definiteness of the estimated between-study (co)variance matrix.
#'  
#' @return A list with components named as the arguments.
#'
#' @examples
#' ## Loading data
#' data("alcohol_cvd")
#' 
#' ## print the iterations (see ?optim) and change the default for starting values
#' dosresmeta(formula = logrr ~ dose, type = type, id = id, se = se, 
#'            cases = cases, n = n, data = alcohol_cvd, proc = "1stage",
#'            control = list(showiter = TRUE, igls.iter = 20))
#'            
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{dosresmeta-package}}, \code{\link{mvmeta.control}}
#' 
#' @references 
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. 
#' Statistics in Medicine, 31(29), 3821-3839.
#' 
#' @export dosresmeta.control

dosresmeta.control <- function(optim = list(), showiter = FALSE, maxiter = 1000, 
                               initPsi = NULL, igls.iter = 10, gr = FALSE,
                               reltol = sqrt(.Machine$double.eps),
                               set.negeigen = sqrt(.Machine$double.eps)){
   optim <- modifyList(list(fnscale = -1, maxit = maxiter, reltol = reltol), 
                       optim)
   if (showiter){
      optim$trace <- 6
      optim$REPORT <- 1
   }
   ## gradient still to be defined (if actually needed)
   list(optim = optim, showiter = showiter, maxiter = maxiter, 
        initPsi = initPsi, igls.iter = igls.iter, gr = gr,
        reltol = reltol, set.negeigen = set.negeigen)
}


#' Fixed-Effects Estimator for dosresmeta Models
#' 
#' @description This function implements a generalized least square estimator for fixed-effects dose-response meta-analysis. It is meant to be used internally 
#' and not directly run by the users.
#' 
#' @param Xlist a m-dimensional list of study-specific design matrices for the fixed-effects part of the model.
#' @param Zlist a m-dimensional list of study-specific design matrices for the random-effects part of the model.
#' @param ylist a m-dimensional list of study-specific of vectors of estimated outcomes.
#' @param Slist a m-dimensional list of within-study (co)variance matrices of estimated outcomes.
#' @param nalist a m-dimensional list of k-dimensional study-specific logical vectors, identifying missing outcomes.
#' @param q,nall numeric scalars: number of predictors, number of observations (excluding missing).
#' @param control list of parameters for controlling the fitting process, usually internally set to default values by \code{dosresmeta.control}.
#' @param \dots further arguments passed to or from other methods. Currently not used.
#' 
#' @details The estimation involves only the \eqn{p} fixed-effects coefficients. The routine is based on a standard generalized least square (GLS) algorithm 
#' implemented in the internal function glsfit. The between-study (co)variance matrix is set to zero, so the marginal (co)variance matrix, composed only by 
#' elements of the within-study component, is assumed as completely known. Similarly to the likelihood-based estimators implemented in 
#' \code{\link{dosresmeta.ml}} and \code{\link{dosresmeta.reml}}, the computation involves Cholesky and and QR decompositions for computational stability and efficiency.
#' 
#' @return This function returns an intermediate list object, whose components are then processed by \code{\link{dosresmeta.fit}}. Other components are 
#' added later through mvmeta to finalize an object of class "\code{dosresmeta}".
#' 
#' @references
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. 
#' Statistics in Medicine, 31(29), 3821-3839.
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{dosresmeta-package}}, \code{\link{dosresmeta.ml}}
#' 
#' @examples
#' data("alcohol_crc")
#' 
#' ## Fixed-effect dose-response model assuming linearity
#' dosresmeta(formula = logrr ~ dose, type = type, id = id, se = se, 
#'            cases = cases, n = peryears, data = alcohol_crc, method = "fixed")
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @name dosresmeta.fixed

dosresmeta.fixed <- function (Xlist, Zlist, ylist, Slist, nalist, q, nall, control, ...){
   Psi <- diag(0, q)
   gls <- glsfit(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   qrinvtUX <- qr(gls$invtUX)
   R <- qr.R(qrinvtUX)
   #Qty <- qr.qty(qrinvtUX, gls$invtUy)
   vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtUX))))
   res <- NULL
   fitted <- lapply(Xlist, "%*%", gls$coef)
   rank <- qrinvtUX$rank
   pconst <- -0.5 * nall * log(2 * pi)
   pres <- -0.5 * (crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
   pdet <- -sum(sapply(gls$Ulist, function(U) sum(log(diag(U)))))
   logLik <- as.numeric(pconst + pdet + pres)
   list(coefficients = gls$coef, vcov = vcov, residuals = res, 
        fitted.values = fitted, df.residual = nall - rank, rank = rank, 
        logLik = logLik, control = control)
}


#' ML and REML Estimators for dosresmeta Models
#' 
#' @description These functions implement maximum likeliihood (ML) and restricted maximum likelihood (REML) estimators for random-effects dose-response 
#' meta-analysis. They are meant to be used internally and not directly run by the users.
#' 
#' @inheritParams dosresmeta.fixed
#'  
#' @details The estimation involves \eqn{p} fixed-effects coefficients and the \eqn{p(p+1)/2} random-effects parameters defining the between-study (co)variance matrix.
#' The hybrid estimation procedure is based first on few runs of iterative generalized least square algorithm and then quasi-Newton iterations, 
#' using specific likelihood functions, until convergence. The estimation algorithm adopts a profiled (or concentrated) approach, that is expressed 
#' only in terms of the random-effects parameters. Cholesky and and QR decompositions are used for computational stability and efficiency, and for assuring the 
#' positive-definiteness of the estimated between-study (co)variance matrix. See the help page for the likelihood functions for further details. 
#' 
#' @return These functions return an intermediate list object, whose components are then processed by \code{\link{dosresmeta.fit}}. Other components are added later 
#' through \code{\link{dosresmeta}} to finalize an object of class "\code{dosresmeta}".
#' 
#' @references
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. 
#' Statistics in Medicine, 31(29), 3821-3839.
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{dosresmeta-package}}, \code{\link{dosresmeta.ml}}
#' 
#' @examples
#' 
#' data("alcohol_cvd")
#' 
#' ## Random-effect dose-response model assuming linearity, ML estimator
#' lin.ml <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                      se = se, cases = cases, n = n, data = alcohol_cvd,
#'                      , method = "ml")
#' summary(lin.ml)
#' 
#' ## Random-effect dose-response model assuming linearity, REML estimator
#' lin.reml <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
#'                        se = se, cases = cases, n = n, data = alcohol_cvd,
#'                        , method = "reml")
#' summary(lin.reml)
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @name dosresmeta.ml

dosresmeta.ml <- function(Xlist, Zlist, ylist, Slist, nalist, q, nall, control, ...){
   par <- initpar(Xlist, Zlist, ylist, Slist, nalist, q, control)
   fn <- mlprof.fn
   gr <- NULL
   opt <- optim(par = par, fn = fn, gr = gr, Xlist = Xlist, Zlist = Zlist, 
                ylist = ylist, Slist = Slist, nalist = nalist, q = q, 
                nall = nall, ctrl = control, 
                method = "Nelder-Mead", control = control$optim)
   Psi <- par2Psi(opt$par, q)
   gls <- glsfit(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   qrinvtUX <- qr(gls$invtUX)
   R <- qr.R(qrinvtUX)
   #    Qty <- qr.qty(qrinvtUX, gls$invtUy)
   vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtUX))))
   res <- NULL
   fitted <- lapply(Xlist, "%*%", gls$coef)
   rank <- qrinvtUX$rank
   list(coefficients = gls$coef, vcov = vcov, Psi = Psi, residuals = res, 
        fitted.values = fitted, df.residual = nall - rank - length(par), 
        rank = rank, logLik = opt$value, converged = opt$convergence == 
           0, niter = opt$counts[[1]], control = control)
}


#' @rdname dosresmeta.ml

dosresmeta.reml <- function(Xlist, Zlist, ylist, Slist, nalist, q, nall, control, ...){
   par <- initpar(Xlist, Zlist, ylist, Slist, nalist, q, control)
   fn <- remlprof.fn
   gr <- NULL
   opt <- optim(par = par, fn = fn, gr = gr, Xlist = Xlist, Zlist = Zlist, 
                ylist = ylist, Slist = Slist, nalist = nalist, q = q, 
                nall = nall, ctrl = control, 
                method = "Nelder-Mead", control = control$optim)
   Psi <- par2Psi(opt$par, q)
   gls <- glsfit(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   qrinvtUX <- qr(gls$invtUX)
   R <- qr.R(qrinvtUX)
   #    Qty <- qr.qty(qrinvtUX, gls$invtUy)
   vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtUX))))
   res <- NULL
   fitted <- lapply(Xlist, "%*%", gls$coef)
   rank <- qrinvtUX$rank
   list(coefficients = gls$coef, vcov = vcov, Psi = Psi, residuals = res, 
        fitted.values = fitted, df.residual = nall - rank - length(par), 
        rank = rank, logLik = opt$value, converged = opt$convergence == 
           0, niter = opt$counts[[1]], control = control)
}


#' Likelihood Functions for dosresmeta Models
#' 
#' @description These functions compute the value of the log-likelihood for random-effects dose-response meta-analysis,
#'  in terms of model parameters. They are meant to be used internally and not directly run by the users.
#' 
#' @param par a vector representing the random-effects parameters defining the between-study (co)variance matrix.
#' @param Psi a p x p matrix representing the current estimate of the between-study (co)variance matrix.
#' @param Xlist a m-dimensional list of study-specific design matrices for the fixed-effects part of the model.
#' @param Zlist a m-dimensional list of study-specific design matrices for the random-effects part of the model.
#' @param ylist a m-dimensional list of study-specific of vectors of estimated outcomes.
#' @param Slist a m-dimensional list of within-study (co)variance matrices of estimated outcomes.
#' @param nalist a m-dimensional list of k-dimensional study-specific logical vectors, identifying missing outcomes.
#' @param p,q,nall numeric scalars: number of predictors, number of observations (excluding missing).
#' @param ctrl list of parameters for controlling the fitting process, usually internally set to default values by 
#' \code{dosresmeta.control}.
#'
#' @details These functions are called internally by the fitting functions \code{\link{dosresmeta.ml}} and \code{\link{dosresmeta.reml}} to 
#' perform iterative optimization algorithms for estimating random effects meta-analytical models.
#' 
#' The maximization of the (restricted) likelihood starts with few runs of an iterative generalized least square algorithm implemented in \code{iter.igls}. 
#' This can be regarded as a fast and stable way to get starting values close to the maximum for the Quasi-Newton iterative algorithm, implemented in 
#' \code{\link{optim}}. Alternatively, starting values can be provided by the user in the control list (see \code{\link{mvmeta.control}}). 
#' 
#' These functions actually specify the profiled version of the (restricted) likelihood, expressed only in terms of random-effects parameters, while the 
#' estimate of the fixed-effects coefficients is provided at each iteration by the internal function \code{glsfit}, based on the current value of 
#' the between-study (co)variance matrix. At convergence, the value of this profiled version is identical to the full (restricted) likelihood. 
#' This approach is computationally efficient, as it reduces the number of parameters in the optimization routine.
#' 
#' The parameterization of the between-study (co)variance matrix ensures the positive-definiteness of the estimated matrix. A Cholesky decomposition is then 
#' performed on the marginal (co)variance matrix in order to re-express the problem as standard least square equations, an approach which speeds up the 
#' computation of matrix inverses and determinants. These equations are finally solved through a QR decomposition, which guarantees stability. 
#' 
#' @return \code{mlprof.fn} and \code{remlprof.fn} return the value of the (restricted) log-likelihood for a given set of 
#' parameters in \code{par}. \code{iter.igls} returns an updated estimate of \code{Psi} given its initial value or the value at 
#' the previous iteration.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{dosresmeta}, \code{\link{mvmeta.fit}}, \code{\link{dosresmeta.control}}, \code{\link{mlprof.fn}}
#' 
#' @name mlprof.fun
remlprof.fn <- function(par, Xlist, Zlist, ylist, Slist, nalist, q, nall, ctrl){ 
   Psi <- par2Psi(par, q)
   gls <- glsfit(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   pconst <- -0.5 * (nall - length(gls$coef)) * log(2 * pi)
   pres <- -0.5 * (crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
   pdet1 <- -sum(sapply(gls$Ulist, function(U) sum(log(diag(U)))))
   tXWXtot <- sumlist(lapply(gls$invtUXlist, crossprod))
   pdet2 <- -sum(log(diag(chol(tXWXtot))))
   as.numeric(pconst + pdet1 + pdet2 + pres)
}

#' @name mlprof.fun
## Not actually used
remlprof.gr <- function(par, Xlist, ylist, Slist, nalist, p, nall, ctrl){
   L <- diag(0, p)
   L[lower.tri(L, diag = TRUE)] <- par
   U <- t(L)
   Psi <- crossprod(U)
   gls <- glsfit(Xlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   tXWXlist <- lapply(gls$invtUXlist, crossprod)
   tXWXtot <- sumlist(lapply(gls$invtUXlist, crossprod))
   invtXWXtot <- chol2inv(chol(tXWXtot))
   invSigmalist <- lapply(gls$invUlist, tcrossprod)
   reslist <- mapply(function(X, y) y - X %*% gls$coef, Xlist, 
                     ylist, SIMPLIFY = FALSE)
   ind1 <- rep(1:p, p:1)
   ind2 <- unlist(sapply(1:p, seq, to = p))
   gradchol.reml(par, U, invtXWXtot, ind1, ind2, Xlist, invSigmalist, 
                 reslist, nalist, p)
}


#' @rdname mlprof.fun
mlprof.fn <- function (par, Xlist, Zlist, ylist, Slist, nalist, q, nall, ctrl){
   Psi <- par2Psi(par, q)
   gls <- glsfit(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   pconst <- -0.5 * nall * log(2 * pi)
   pres <- -0.5 * (crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
   pdet <- -sum(sapply(gls$Ulist, function(U) sum(log(diag(U)))))
   as.numeric(pconst + pdet + pres)
}

#' @rdname mlprof.fun
## Not actually used
mlprof.gr <- function(par, Xlist, ylist, Slist, nalist, p, nall, ctrl){
   L <- diag(0, p)
   L[lower.tri(L, diag = TRUE)] <- par
   U <- t(L)
   Psi <- crossprod(U)
   gls <- glsfit(Xlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   invSigmalist <- lapply(gls$invUlist, tcrossprod)
   reslist <- mapply(function(X, y) as.numeric(y - X %*% gls$coef), 
                     Xlist, ylist, SIMPLIFY = FALSE)
   ind1 <- rep(1:p, p:1)
   ind2 <- unlist(sapply(1:p, seq, to = p))
   gradchol.ml(par, U, ind1, ind2, invSigmalist, reslist, Xlist, nalist, p)
}

#' @rdname mlprof.fun
iter.igls <- function(Psi, Xlist, Zlist, ylist, Slist, nalist, q){
   gls <- glsfit(Xlist, Zlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
   npar <- q * (q + 1)/2
   #indMat <- xpndMat(seq(npar))
   eSigmalist <- lapply(gls$Sigmalist, function(Sigma) Sigma %x% 
                           Sigma)
   flist <- mapply(function(y, S, X) {
      return(as.numeric(tcrossprod(y - X %*% gls$coef)) - as.numeric(S))
   }, ylist, Slist, Xlist, SIMPLIFY = FALSE)
   hlist <- lapply(Xlist, function(y){
      lapply(seq(npar), function(x){
         m <- matrix(0, q, q)
         #m[!upper.tri(m)][x] <- 1
         m[!upper.tri(m)][x] <- m[!lower.tri(m)][x] <- 1
         m
      })
   })
   Zylist <- mapply(function(x, h){
      ris <- list()
      for (i in 1:length(h)){
         ris <- c(ris, list(as.numeric(x %*% h[[i]] %*% t(x))))
      }
      return(do.call("cbind", ris))
   }, Zlist, hlist, SIMPLIFY = FALSE)
   eUlist <- lapply(eSigmalist, chol)
   inveUlist <- lapply(eUlist, function(U) backsolve(U, diag(ncol(U))))
   invteUZylist <- mapply(function(inveU, Zy) 
      crossprod(inveU, Zy), inveUlist, Zylist, SIMPLIFY = FALSE)
   invteUflist <- mapply(function(inveU, f) 
      crossprod(inveU, f), inveUlist, flist, SIMPLIFY = FALSE)
   invteUZy <- do.call("rbind", invteUZylist)
   invteUf <- do.call("rbind", invteUflist)
   #theta <- unique(as.numeric(qr.solve(invteUZ, invteUf)))
   theta <- as.numeric(qr.solve(invteUZy, invteUf))
   Psi <- xpndMat(theta)
   eig <- eigen(Psi)
   eig$values <- pmax(eig$values, 10^-8)
   eig$vectors %*% diag(eig$values, q) %*% t(eig$vectors)
}