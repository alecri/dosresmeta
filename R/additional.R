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
