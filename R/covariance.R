#' Computes mean and standardized mean differences for continuous outcome with corresponding
#' co(variance) matrix
#' 
#' @description This internal function computes mean and standardized mean of a continuous outcome with the corresponding 
#' variances. It also reconstructs the covariance matrix from the available data.
#' 
#' @param y a vector defining the mean outcome for each treatment level.
#' @param sd a vector defining the standard deviation of the outcome for each treatment level.
#' @param n a vector defining the number of subjects for each treatment level.
#' @param measure character string, indicating the measure to be calculated. Options are \code{md}
#' and \code{smd} for mean difference and standardized mean difference, respectively.
#' @param method character string indicating the method to be used. Options are \code{cohens}, \code{hedges}, and \code{glass}. 
#' @param data an optional data frame (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the previous arguments.
#' 
#' @return A list containing the following
#' \tabular{ll}{
#' \code{y} \tab mean or standardized mean differences for each treatment level,
#' included the referent one (0 by calculation).\cr
#' \code{v} \tab variances corresponding to the mean or standardized mean differences for each 
#' treatment level, included the referent one (0 by calculation)\cr
#' \code{S} \tab co(variance) matrix for the non-referent mean or standardized mean differences.\cr
#' }
#' 
#' @details This is an internal function called by \code{\link{dosresmeta}} to reconstruct the (co)variance matrix of the 
#' outcome variable. The function is expected to be extended and/or modified at every release of the package
#'
#' @seealso \code{\link{covar.logrr}}, \code{\link{dosresmeta}} 
#'
#' @examples
#' ## Loading the data
#' data("ari")
#' 
#' ## Obtaining standardized mean differences, variances, and (co)varinace 
#' ## matrix for the first study (id = 1)
#' covar.smd(y, sd, n, measure = "smd", data = subset(ari, id == 1))
#' 
#' ## Obtaining mean differences, variances, and (co)varinace matrices for the all the studies
#' cov.md <- by(ari, ari$id, function(x) covar.smd(y, sd, n, "md", data = x))
#' 
#' ## Extracting mean differences
#' unlist(lapply(cov.md, function(x) x$y))
#' ## Extracting variances for the mean differences
#' unlist(lapply(cov.md, function(x) x$v))
#' ## List of the (co)variance matrices for the mean differences
#' lapply(cov.md, function(x) x$S)
#'  
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references 
#' 
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2009). The handbook of 
#' research synthesis and meta-analysis. Russell Sage Foundation.
#' 
#' @export covar.smd

covar.smd <- function(y, sd, n, measure = "md", method = "cohens", data){
   if (missing(data)) 
      data <- NULL
   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   }
   else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }
   method <- match.arg(method, c("cohens", "hedges", "glass"))
   mf <- match.call()
   mf.y <- mf[[match("y", names(mf))]]
   mf.sd <- mf[[match("sd", names(mf))]]
   mf.n <- mf[[match("n", names(mf))]]
   y <- eval(mf.y, data, enclos = sys.frame(sys.parent()))
   sd <- eval(mf.sd, data, enclos = sys.frame(sys.parent()))
   n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
   y <- scale(y, y[1], F)
   sdPooled <- sqrt(sum((n-1)*sd^2)/sum(n-1))
   if (measure == "md"){
      v <- sdPooled^2 * (n + n[1]) / (n * n[1])
      v[1] <- 0
      S <- matrix(sd[1]^2/n[1], ncol = length(v[v != 0]), 
                  nrow = length(v[v != 0]))
      diag(S) <- v[v != 0]
   }
   if (measure == "smd"){
      if (method == "cohens"){
         y <- y/sdPooled
         S <- 1/n[1] + tcrossprod(y[-1])/(2*sum(n))
      }
      if (method == "hedges"){
         y <- (y/sdPooled) * (1 - (3/((4 * sum(n)) - 9)))
         S <- 1/n[1] + tcrossprod(y[-1])/(2 * (sum(n) - 3.94))
      }
      else if (method == "glass"){
         y <- y/sd[1]
         S <- 1/n[1] + tcrossprod(y[-1])/(2 * (n[1] - 1))
      }
      v <- c(0, 1/n[-1] + diag(S))
      diag(S) <- v[-1]
   } 
   list(y = y, v = v, S = S)
}


#' Computes the covariance matrix for a set of log relative risks
#' 
#' @description Reconstructs the covariance matrix for a set of (reported) log relative risks, given the number of cases and 
#' the number of total persons or person-years for each treatment (dose) level.
#' 
#' @inheritParams hamling
#' @param covariance method to approximate the coviariance among set of reported log relative risks, "\code{gl}" for the method proposed by Greenland and Longnecker  
#' (default), "\code{h}" for the method proposed by Hamling.
#'  
#' @return The (co)variance matrix of the log relative risks.
#' 
#' @details This is an internal function called by \code{\link{dosresmeta}} to reconstruct the (co)variance matrix of the (adjusted) log relative risks. The function
#' calls, depending on the choosen method, \code{\link{grl}} (default) or \code{\link{hamling}} to reconstruct the effective counts corresponding to the (adjusted) log
#' relative risks as well as their standard errors. From these it computes the covariance matrix; analytical formulas can be found in the referenced article.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: 
#' examples, an evaluation of approximations, and software. American journal of epidemiology, 175(1), 66-73.
#' 
#' @seealso \code{\link{grl}}, \code{\link{hamling}}, \code{\link{covar.smd}}, \code{\link{dosresmeta}}
#' 
#' @examples
#' ## Loading data
#' data("alcohol_cvd")
#' 
#' ## Obtaining the (co)variance matrix of log RR for the first study (id = 1)
#' covar.logrr(y = logrr, v = I(se^2), cases = cases, n = n, type = type, 
#'             data = subset(alcohol_cvd, id == 1))
#' 
#' ## Obtaining the (co)variance matrices of log RRfor all study
#' by(alcohol_cvd, alcohol_cvd$id, function(x)
#'    covar.logrr(y = logrr, v = I(se^2), cases = cases, n = n, 
#'                type = type, data = x))
#' 
#' ## Restructuring the previous results in a list of matrices
#' do.call("list", by(alcohol_cvd, alcohol_cvd$id, function(x)
#'    covar.logrr(y = logrr, v = I(se^2), cases = cases, n = n, type = type,
#'                data = x)))
#' 
#'@export covar.logrr

covar.logrr <- function(cases, n, y, v, type, data, covariance = "gl"){
   if (missing(data)) 
      data <- NULL
   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   }
   else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }
   covariance <- match.arg(covariance, c("gl", "hamling"))
   mf <- match.call(expand.dots = FALSE)
   mf.cases <- mf[[match("cases", names(mf))]]
   mf.n <- mf[[match("n", names(mf))]]
   mf.y <- mf[[match("y", names(mf))]]
   mf.v <- mf[[match("v", names(mf))]]
   mf.type <- mf[[match("type", names(mf))]]
   cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
   n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
   y <- eval(mf.y, data, enclos = sys.frame(sys.parent()))
   v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
   v[is.na(v)] <- 0
   type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
   if (is.null(type)) 
      type <- as.vector(mf.type)
   ps <- if (covariance == "gl")
      grl(y, v, cases, n, type)
   else
      hamling(y, v, cases, n, type)
   rcorr <- switch(as.character(type[1]),
                   cc = {
                      s0 <- 1/ps[v==0, 1] + 1/(ps[v==0, 2] - ps[v==0, 1])
                      si <- s0 + 1/ps[v!=0, 1] + 1/(ps[v!=0, 2] - ps[v!=0, 1])
                      s0/tcrossprod(si)^0.5
                   },
                   ir = {
                      s0 <- 1/ps[v==0, 1]
                      si <- s0 + 1/ps[v!=0, 1]
                      s0/tcrossprod(si)^0.5
                   },
                   ci = {
                      s0 <- 1/ps[v==0, 1] - 1/ps[v==0, 2]
                      si <- s0 + 1/ps[v!=0, 1] - 1/ps[v!=0, 2]
                      s0/tcrossprod(si)^0.5
                   })
   diag(rcorr) <- 1
   tcrossprod(v[v != 0])^0.5 * rcorr
}


#' Approximating effective-counts as proposed by Hamling
#' 
#' @description Reconstructs the set of pseudo-numbers (or "effective" numbers) of cases and non-cases consistent
#' with the input data (log relative risks). The method was first proposed in 2008 by Hamling.
#' 
#' @param y a vector, defining the (reported) log relative risks.
#' @param v a vector, defining the variances of the reported log relative risks.
#' @param cases a vector, defining the number of cases for each exposure level.
#' @param n a vector, defining the total number of subjects for each exposure level. For incidence-rate data \code{n} indicates the amount of person-time within 
#' each exposure level.
#' @param type a vector (or a character string), specifying the design of the study. Options are
#' \code{cc}, \code{ir}, and \code{ci}, for case-control, incidence-rate, and cumulative incidence data, respectively.
#' @param data an optional data frame (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the previous arguments.
#' 
#' @return A list containing the following
#' \tabular{ll}{
#' \code{y} \tab mean or standardized mean differences for each treatment level,
#' included the referent one (0 by calculation).\cr
#' \code{v} \tab variances corresponding to the mean or standardized mean differences for each 
#' treatment level, included the referent one (0 by calculation)\cr
#' \code{S} \tab co(variance) matrix for the non-referent mean or standardized mean differences.\cr
#' } 
#' 
#' @details The function reconstructs the effective counts corresponding to the multivariable adjusted log relative risks as well as their standard errors. 
#' A unique solution is guaranteed by keeping the ratio non-cases to cases and the fraction of unexposed subjects equal to the unadjusted data (Hamling). 
#' See the referenced article for a complete description of the algorithm implementation.
#'
#' @examples
#' ## Loading data
#' data("alcohol_cvd")
#' 
#' ## Obtaining pseudo-counts for the first study (id = 1)
#' hamling(y = logrr, v = I(se^2), cases = cases, n = n, type = type, 
#' data = subset(alcohol_cvd, id == 1))
#'    
#' ## Obtaining pseudo-counts for all study
#' by(alcohol_cvd, alcohol_cvd$id, function(x)
#' hamling(y = logrr, v = I(se^2), cases = cases, n = n, type = type, data = x))
#'
#' ## Restructuring the previous results in a matrix
#' do.call("rbind", by(alcohol_cvd, alcohol_cvd$id, function(x)
#'    hamling(y = logrr, v = I(se^2), cases = cases, n = n, type = type,
#'       data = x)))
#'  
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{\link{covar.logrr}}, \code{\link{grl}}
#' 
#' @references 
#' Hamling, J., Lee, P., Weitkunat, R., Ambuhl, M. (2008). Facilitating meta-analyses by deriving relative effect and precision estimates for alternative 
#' comparisons from a set of estimates presented by exposure level or disease category. Statistics in medicine, 27(7), 954-970.
#' 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, an evaluation 
#' of approximations, and software. American journal of epidemiology, 175(1), 66-73.
#' 
#' @export hamling

hamling <- function(y, v, cases, n, type, data){
   if (missing(data)) 
      data <- NULL
   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   }
   else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }
   mf <- match.call(expand.dots = FALSE)
   mf.y <- mf[[match("y", names(mf))]]
   mf.v <- mf[[match("v", names(mf))]]
   mf.cases <- mf[[match("cases", names(mf))]]
   mf.n <- mf[[match("n", names(mf))]]
   mf.type <- mf[[match("type", names(mf))]]
   y <- eval(mf.y, data, enclos = sys.frame(sys.parent()))
   v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
   v[is.na(v)] <- 0
   cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
   n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
   type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
   if (is.null(type)) 
      type <- as.vector(mf.type)
   ## consistent with glst2.ado and ecov.ado
   p0 <- if (as.character(type[1]) == "cc")
      (n - cases)[v==0]/sum(n - cases)
   else
      n[v==0]/sum(n)
   z0 <- if (as.character(type[1]) == "cc")
      sum(n - cases)/sum(cases)
   else
      sum(n)/sum(cases)
   init <- c(cases[v==0], n[v==0])
   opt <- optim(init, fun.h, v = v, y = y, type = type, p0 = p0, z0 = z0)
   pscounts <- est.ps.h(opt$par, v, y, type)
   pscounts
}

#' @noRd
est.ps.h <- function(param, v, y, type){
   A0 <- param[1]
   N0 <- param[2]
   ps <- switch(as.character(type[1]), 
                cc = cbind(A <- (1+exp(y)*A0/(N0 - A0))/(v-1/A0-1/(N0 - A0)),
                           N = A + (1+(N0 - A0)/(A0*exp(y)))/(v-1/A0-1/(N0 - A0))),
                ir = cbind(A = (1 - exp(y)*A0/N0)/(v - 1/A0 + 1/N0),
                           N = (N0/(A0*exp(y))-1)/(v-1/A0 + 1/N0)),
                ci = cbind(A = 1/(v - 1/A0),
                           N = (N0/(A0*exp(y)))/(v-1/A0)))
   ps[v==0, ] <- param
   colnames(ps) <- c("A", "N") 
   return(ps)
}

#' @noRd
fun.h <- function(par, v, y, type, p0, z0){
   ps <- est.ps.h(par, v, y, type)
   if (as.character(type[1]) == "cc"){
      p1 <- (ps[v==0, 2] - ps[v==0, 1])/(sum(ps[, 2] - ps[, 1]))
      z1 <- sum(ps[, 2] - ps[, 1])/sum(ps[, 1])
   }
   else {
      p1 <- ps[v==0, 2]/sum(ps[, 2])
      z1 <- sum(ps[, 2])/sum(ps[, 1])
   }
   ((p1 - p0)/p0)^2+((z1 - z0)/z0)^2
}

#' Approximating effective-counts as proposed by Greenland & Longnecker
#' 
#' @description Reconstructs the set of pseudo-numbers (or 'effective' numbers) of cases and non-cases consistent
#' with the input data (log relative risks). The method was first proposed in 1992 by Greenland and Longnecker.
#'
#' @inheritParams hamling 
#' @param tol define the tolerance.
#' 
#' @details The function reconstructs the effective counts corresponding to the multivariable adjusted log relative risks as well as their standard errors. 
#' A unique solution is guaranteed by keeping the margins of the table of pseudo-counts equal to the margins of the crude or unadjusted data 
#' (Greenland and Longnecker 1992). See the referenced article for a complete description of the algorithm implementation.
#' 
#' @return The results are returned structured in a matrix
#' \tabular{ll}{
#' \code{A} \tab approximated number of effective cases. \cr
#' \code{N} \tab approximated total number of effective subjects. \cr
#' } 
#' 
#' @examples
#' ## Loading data
#' data("alcohol_cvd")
#' 
#' ## Obtaining pseudo-counts for the first study (id = 1)
#' grl(y = logrr, v = I(se^2), cases = cases, n = n, type = type,
#'    data = subset(alcohol_cvd, id == 1))
#'    
#' ## Obtaining pseudo-counts for all study
#' by(alcohol_cvd, alcohol_cvd$id, function(x)
#'    grl(y = logrr, v = I(se^2), cases = cases, n = n, type = type, data = x))
#'
#' ## Restructuring the previous results in a matrix
#' do.call("rbind", by(alcohol_cvd, alcohol_cvd$id, function(x)
#'    grl(y = logrr, v = I(se^2), cases = cases, n = n, type = type, data = x)))
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @seealso \code{\link{covar.logrr}}, \code{\link{hamling}}
#' 
#' @references 
#' Greenland, S., Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. American journal of epidemiology, 135(11), 1301-1309.
#' 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, an evaluation of approximations, and software. 
#' American journal of epidemiology, 175(1), 66-73.
#' 
#' @export grl

grl <- function(y, v, cases, n, type, data, tol = 1e-05){
   if (missing(data)) 
      data <- NULL
   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   }
   else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }
   mf <- match.call(expand.dots = FALSE)
   mf.y <- mf[[match("y", names(mf))]]
   mf.v <- mf[[match("v", names(mf))]]
   mf.cases <- mf[[match("cases", names(mf))]]
   mf.n <- mf[[match("n", names(mf))]]
   mf.type <- mf[[match("type", names(mf))]]
   y <- eval(mf.y, data, enclos = sys.frame(sys.parent()))
   v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
   v[is.na(v)] <- 0
   cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
   n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
   type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
   if (is.null(type)) 
      type <- as.vector(mf.type)
   Ax <- Axp <- cases
   repeat{
      A0 <- sum(cases) - sum(Ax[v!=0])
      cx <- if (!as.character(type[1]) == "ir")
         1/Ax + 1/(n - Ax) 
      else
         1/Ax
      e <- if (!as.character(type[1]) == "ir")
         y[v!=0] + log(A0) + log(n[v!=0]-Ax[v!=0]) - log(Ax[v!=0]) - 
         log(n[v==0] - A0) 
      else 
         y[v!=0] + log(A0) + log(n[v!=0]) - log(Ax[v!=0]) - log(n[v==0])
      H <- diag(cx[v!=0] + cx[v==0], nrow = sum(v!=0))
      H[upper.tri(H)] <- H[lower.tri(H)] <- cx[v==0]
      Axp[v==0] <- A0
      Axp[v!=0] <- Ax[v!=0] + solve(H) %*% e
      delta <- sum((Axp[v!=0] - Ax[v!=0])^2)
      if (delta < tol)
         break
      Ax <- Axp
   }
   cbind(A = Axp, N = n)
}

#' @noRd
change_ref <- function(y, v, cases, n, type, data, ref = 1,
                       method = "hamling", expo = FALSE){
   if (missing(data)) 
      data <- NULL
   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   }
   else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }
   mf <- match.call(expand.dots = FALSE)
   mf.y <- mf[[match("y", names(mf))]]
   mf.v <- mf[[match("v", names(mf))]]
   mf.cases <- mf[[match("cases", names(mf))]]
   mf.n <- mf[[match("n", names(mf))]]
   mf.type <- mf[[match("type", names(mf))]]
   y <- eval(mf.y, data, enclos = sys.frame(sys.parent()))
   v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
   v[is.na(v)] <- 0
   cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
   n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
   type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
   if (is.null(type)) 
      type <- as.vector(mf.type)
   ps <- if (method == "hamling"){
      hamling(y, v, cases, n, type)
   } else if (method == "gl"){
      grl(y, v, cases, n, type)
   }
   rr <- if (type[1] == "cc"){
      (ps[, 1]*(ps[ref, 2] - ps[ref, 1]))/(ps[ref, 1]*(ps[, 2] - ps[, 1]))
   } else {
      (ps[, 1]/ps[, 2])/(ps[ref, 1]/ps[ref, 2])
   }
   var_logrr <- if (type[1] == "cc"){
      1/ps[, 1] + 1/(ps[, 2] - ps[, 1]) + 1/ps[ref, 1] + 1/(ps[ref, 2] - ps[ref, 1])
   } else if (type[1] == "ci") {
      1/ps[, 1] - 1/ps[, 2] + 1/ps[ref, 1] - 1/ps[ref, 2]
   } else if (type[1] == "ir") {
      1/ps[, 1] + 1/ps[, 2]
   }
   var_logrr[ref] <- 0
   lb_rr <- exp(log(rr) - qnorm(.975)*sqrt(var_logrr))
   ub_rr <- exp(log(rr) + qnorm(.975)*sqrt(var_logrr))
   dat <- if (expo){
      data.frame(ps, rr = rr, lb_rr = lb_rr, ub_rr = ub_rr)
   } else {
      data.frame(ps, logrr = log(rr), v = var_logrr, 
                 loglbrr = log(lb_rr), logubrr = log(ub_rr))
   }
   colnames(dat) <- paste(colnames(dat), ref, sep = ".")
   return(dat)
}