.onAttach <- function(libname, pkgname) {
   packageStartupMessage("This is dosresmeta 2.0.0. For an overview type: help('dosresmeta-package').")
}

#' Performing multivariate dose-response meta-analysis  
#' @name dosresmeta-package
#' @docType package
#' @import mvmeta
#' 
#' @description The package \code{dosresmeta} consists of a collection of functions to estimate dose-response 
#' relations from summarized dose-response data for both continuous and binary outcomes, and to combine them according to principles of (multivariate) random-effects model. 
#' 
#' @section Modeling framework:
#' 
#' Dose-response meta-analysis represents a specific type of meta-analysis. Aim of such analysis is to reconstruct and combine study-specific curves from summarized 
#' dose-response data. Greenland and Longnecker originally developed the methodology in 1992 for pooling associations from epidemiological studies of binary outcomes. 
#' Extensions are currently proposed for other types of outcomes (e.g. continuous) from others study design, such as clinical trials.
#' 
#' The summarized dose-response data are most often presented in a tabular way, reporting the levels of the exposure (doses) and the corresponding outcome variable. 
#' The latter is usually expressed as contrast to the unexposed or baseline category (referent level). Examples are (log) relative risks, (log) odds ratios, 
#' (log) incidence rate ratios, mean differences, and standardized mean differences. Thus the outcome cannot be regarded as independent and a (co)variance matrix needs 
#' to be provided or approximated from the available data.See \code{\link{covar.smd}} and \code{\link{covar.logrr}} for more details.
#' 
#' @section Estimation procedure:
#' 
#' The pooled dose-response association can be estimated using two different approaches. The former consists of a two-stage procedure, where the study-specific trend 
#' are first estimated and then pooled across studies. Assuming \eqn{y_j}{yj} is the vector of non-referent outcome values in each of \eqn{i = 1, \dots, m} studies, and 
#' \eqn{X_i}{Xi} the related matrix of \eqn{p} transformations of the exposure (typically \eqn{p = 1, 2}), the dose-response model can be written as 
#' \deqn{y_i = X_i\beta_i + \epsilon_i}{yi = Xi\betai + \epsiloni} with \eqn{S_i}{Si} = (co)variance of \eqn{\epsilon_i}{\epsiloni} known (available or 
#' reconstructed from the available data). 
#' The \eqn{\beta_i}{\betai} are then combined according to principles of (multivariate) random-effects meta-analytical models 
#' \deqn{\beta_i ~ N ( \beta, V_i + \Psi )}{\betai ~ N( \beta, Vi + \Psi )}
#' where \eqn{V_i}{Vi} and \eqn{\Psi} indicate, respectively, the within study (co)variance (obtained in the first stage analysis) and the between study (co)variance.
#' 
#' The alternative approach, instead, consists of a one-stage (also known as pool-first) procedure. The data are pooled by concatenating the vector \eqn{y_i}{yi} and 
#' vectors (or matrices) \eqn{X_i}{Xi}. The (multivariate) random effects-model can be written as \deqn{y_i = X_i\beta + Z_i\eta_i + \epsilon_i}{yi = X_i\betai + Z_i\etai + \epsiloni} 
#' where \eqn{\beta} represents the fixed-effects parameter, \eqn{\eta_i} the vector (or matrix) of unobserved random-effects for thei-th study, and \eqn{Z_i}{Zi} coincides with \eqn{X_i}{Xi}. 
#' The marginal model has a co(variance) matrix equal to \eqn{\Sigma + Z_i\Psi Z_i^top}{\Sigma + Zi\PsiZi'}, where \eqn{\Sigma} is the block diagonal (co)variance with i-th diagonal block \eqn{S_i}{Si}.
#' 
#' The two approaches provide similar results, despite the two-stage procedure may be more stable and faster in terms of convergence. In both the procedures the 
#' aim is to estimate the coefficients \eqn{\beta} and, for random-effects models, the components of the between-study (co)variance matrix \eqn{Psi}. 
#' Different estimators are implemented in the package. The estimation options available are
#' \itemize{
#' \item Fixed-effects
#' \item Maximum likelihood (ML)
#' \item Restricted maximum likelihood (REML)
#' \item Method of moments (currently available only for the two-stage procedure)
#' }
#' 
#' The fixed-effects model is fitted through generalized least squares (GLS), assuming the (co)variance structure, composed by the within-study error only, 
#' as completely known. Among random-effects models, ML and REML approaches provides fit criteria and inferential test derived from likelihood theory, such 
#' as AIC and likelihood ratio test, purticularly useful in a one-stage procedure. Further details on estimation methods are given in the related help pages.
#' 
#' @section Functions and data included in the package:
#'
#' The structure of the package and the internal functions resemble those of the \code{\link{mvmeta}} package. See \code{\link{mvmeta-package}} for a general overview. 
#' The main function is \code{\link{dosresmeta}}, which performs the various models illustrated above. The function returns a list object of class 
#' "\code{dosresmeta}" (see \code{\link{dosresmetaObject}}).
#'
#' The estimation is carried out internally through \code{\link{dosresmeta.fit}}, a wrapper which prepares the data and calls specific estimation functions 
#' for fitting the models, depending on the chosen procedure. For the two-stage procedure, the second part of the analysis is performed using the function \code{\link{mvmeta.fit}}
#' while estimators for random-effects models are implemented in the functions \code{\link{dosresmeta.ml}} and \code{\link{dosresmeta.reml}} for 
#' (restricted) maximum likelihood. For likelihood-based methods, iterative optimizations algorithms are used for maximizing 
#' the (restricted) likelihood. Fitting parameter options are set by \code{\link{dosresmeta.control}}.
#'
#' Method functions are available for objects of class "\code{dosremeta}" (see \code{\link{dosresmetaObject}} for a complete list). The method \code{\link{summary}} 
#' produces a list of class "\code{summary.dosremeta}" for summarizing the fit of the model and providing additional results. The method function \code{\link{predict}} 
#' computes predicted values, optionally for a set of new values of the predictors. \code{\link{blup}} gives the (empirical) best linear unbiased predictions for the unobserved random-effects.
#' Other default or specific method functions for regression can be used on objects of class "\code{dosremeta}", such as \code{logLik}, \code{AIC} and \code{BIC}, among others. 
#' The method function \code{qtest.dosresmeta} (producing an object with class of the same name) performs the Cochran Q test for (residual) heterogeneity currently appropriate only for the two-stage approach.
#'
#' Printing functions for the objects of classes defined above are also provided. Other functions are used internally in the source code, and not exported in the 
#' namespace. For users interested in getting into details of the package structure, these functions can be displayed using the triple colon (':::') operator. 
#' For instance, dosresmeta:::glsfit displays the code of the function glsfit. 
#'
#' The package includes the datasets \code{\link{alcohol_crc}}, \code{\link{alcohol_cvd}}, \code{\link{ari}}, and \code{\link{cc_ex}} as data frames, 
#' which are used in the examples. 
#'   
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' @references
#' Greenland, S., Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. 
#' American journal of epidemiology, 135(11), 1301-1309.
#' 
#' Orsini, N., Bellocco, R.,  Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.
#' 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, 
#' an evaluation of approximations, and software. American journal of epidemiology, 175(1), 66-73.
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. 
#' Statistics in Medicine, 31(29), 3821-3839.
#'   
#' @seealso \code{\link{dosresmeta}} \code{\link{mvmeta}} 

NULL

#' Case-control data on alcohol and breast cancer risk
#'
#' @name cc_ex
#' @description The dataset reports the summarized dose-response results from a case-control study
#' on alcohol and breast cancer, first presented by Rohan and McMichael.
#' @docType data
#' @format A data frame with 4 observations on the following 10 variables:
#' \tabular{ll}{
#' \code{gday} \tab label for exposure levels.\cr
#' \code{dose} \tab assigned dose level.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{control} \tab number of controls for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{crudeor} \tab unadjusted odds ratios for each exposure level.\cr
#' \code{adjrr} \tab adjusted odds ratios for each exposure level.\cr
#' \code{lb} \tab lower bound for the confidence limit of the adjusted odds ratios.\cr
#' \code{ub} \tab upper bound for the confidence limit of the adjusted odds ratios.\cr
#' \code{logrr} \tab natural logarithm of adjusted odds ratios.\cr
#' } 
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references 
#' Rohan, T. E., McMichael, A. J. (1988). Alcohol consumption and risk op breast cancer. International journal of cancer, 41(5), 695-699.
#' 
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. 
#' American journal of epidemiology, 135(11), 1301-1309.
#' 
#' @keywords data
NULL


#' Six published studies on the relation between alcohol intake and vascular disease risk.
#' 
#' @name alcohol_cvd
#' @description The dataset reports the summarized dose-response results from six observational
#' studies on the relation between alcohol intake and vascular disease risk (Qin Liu 2009). Four are case-control studies, 
#' two prospective (cumulative-incidence data). 
#' @docType data
#' @format A data frame with 25 observations on the following 8 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author of the study.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose level.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of adjusted "relative risks".\cr
#' \code{se} \tab standard errornatural for the logarithm of adjusted "relative risks".\cr
#' }
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references 
#' Liu, Q., Cook, N. R., Bergstrom, A., Hsieh, C. C. (2009). A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear 
#' dose-response data. Computational Statistics & Data Analysis, 53(12), 4157-4167.
#'   
#' @keywords data
NULL


#' Eight published studies on the relation between alcohol intake and colon-rectal cancer.
#' 
#' @name alcohol_crc
#' @description The dataset reports the summarized dose-response results from eight prospective
#' #' studies on the relation between alcohol intake and colon-rectal risk (Orsini 2012).
#' @docType data
#' @format A data frame with 48 observations on the following 7 variables:
#' \tabular{ll}{
#' \code{id} \tab label for author's names (id variable).\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose level.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{peryears} \tab amount of person-time for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of adjusted "relative risks".\cr
#' \code{se} \tab standard errornatural for the logarithm of adjusted "relative risks".\cr
#' }
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, 
#' an evaluation of approximations, and software. American journal of epidemiology, 175(1), 66-73.
#'   
#' @keywords data
NULL


#' Five clinical trials on the relation between aripiprazole and schizophrenia
#' 
#' @name ari
#' @description The dataset reports the summarized dose-response results from five clinical trials on the relation between different levels of aripiprazole 
#' and severety of schizophrenia measured usign the PANSS medical score.
#' 
#' @docType data
#' @format A data frame with 18 observations on the following 6 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author of the studies.\cr
#' \code{dose} \tab assigned dose level of aripiprazole (0 for placebo group).\cr
#' \code{y} \tab outcome variable: change in PANNS score after and before treatment.\cr
#' \code{sd} \tab standard deviation of y for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' }
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' @keywords data
NULL