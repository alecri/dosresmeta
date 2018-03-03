#' Case-control data on alcohol and breast cancer risk
#'
#' @name cc_ex
#' @description The dataset reports the summarized dose-response results from a case-control study
#' on alcohol and breast cancer, first presented by Rohan and McMichael.
#' @docType data
#' @format A data frame with 4 observations on the following 10 variables:
#' \tabular{ll}{
#' \code{gday} \tab label for exposure levels.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{control} \tab number of controls for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{crudeor} \tab unadjusted odds ratios for each exposure level.\cr
#' \code{adjrr} \tab adjusted odds ratios for each exposure level.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted odds ratios.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted odds ratios.\cr
#' \code{logrr} \tab natural logarithm of the adjusted odds ratios.\cr
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


#' Cumulative incidence data on high-fat dairy food and colorectal cancer risk
#'
#' @name ci_ex
#' @description The dataset reports the summarized dose-response results from a cumlative-incidence study
#' on high-fat dairy food intake and risk of colorectal cancer, first presented by Larsson, Bergkvist, and Wolk (2005).
#' @docType data
#' @format A data frame with 5 observations on the following 8 variables:
#' \tabular{ll}{
#' \code{dose} \tab assigned dose levels.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{adjrr} \tab adjusted risk ratios for each exposure level.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted risk ratios.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted risk ratios.\cr
#' \code{logrr} \tab natural logarithm of adjusted risk ratios.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk ratios.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#'
#' @references
#' 
#' Larsson, S. C., L. Bergkvist, and A. Wolk. (2005). High-fat dairy food and conjugated 
#' linoleic acid intakes in relation to colorectal cancer incidence in the Swedish Mammography 
#' Cohort. American Journal of Clinical Nutrition 82: 894-900.
#'
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis.
#' American journal of epidemiology, 135(11), 1301-1309.
#'
#' @keywords data
NULL


#' Incidence-rate data on fiber intake and coronary heart disease risk
#'
#' @name ir_ex
#' @description The dataset reports the summarized dose-response results from incidence-rate data 
#' investigating the association between the long-term intake of dietary fiber and 
#' risk of coronary heart disease among women, first presented by Wolk et al. (1999)
#' @docType data
#' @format A data frame with 5 observations on the following 8 variables:
#' \tabular{ll}{
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{adjrr} \tab adjusted incidence rate ratios for each exposure level.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted incidence rate ratios.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted incidence rate ratios.\cr
#' \code{logrr} \tab natural logarithm of adjusted incidence rate ratios.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted incidence rate ratios.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#'
#' @references
#' Wolk, A., J. E. Manson, M. J. Stampfer, G. A. Colditz, F. Hu, F. E. Speizer, C. H. Hennekens, and W. C. Willett. 1999. 
#' Long-term intake of dietary fiber and decreased risk of coronary heart disease among women. Journal of the American Medical Association 281: 1998-2004.
#'
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis.
#' American journal of epidemiology, 135(11), 1301-1309.
#'
#' @keywords data
NULL


#' Nine studies on the relation between milk consumption and ovarian cancer
#'
#' @name milk_ov
#' @description The dataset reports the summarized dose-response results from nine studies on the relation between milk consumption
#' and ovarian cancer.
#'
#' @docType data
#' @format A data frame with 37 observations on the following 12 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ir" or "cc") or person-years (type = "ir") for each exposure level.\cr
#' \code{rr} \tab adjusted risk estimates.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Larsson, S. C., N. Orsini, and A. Wolk. 2005. Milk, milk products and lactose intake and ovarian cancer risk: A meta-analysis of epidemiological studies. 
#' International Journal of Cancer 118: 431-441.
#' 
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis.
#' American journal of epidemiology, 135(11), 1301-1309.
#' 
#' @keywords data
NULL


#' Twenty-two case-control studies on the relation between oral contraceptives use and breast cancer
#'
#' @name oc_breast
#' @description The dataset reports the summarized dose-response results from twenty-two case-control 
#' studies on the relation between oral contraceptives use and breast cancer
#'
#' @docType data
#' @format A data frame with 113 observations on the following 14 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{duration} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ir" or "cc") or person-years (type = "ir") for each exposure level.\cr
#' \code{or} \tab adjusted odds ratios.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted odds ratios.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted odds ratios.\cr
#' \code{logor} \tab natural logarithm of the adjusted odds ratios.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted odds ratios.\cr
#' \code{menopause} \tab indicator variable for a study that included postmenopausal women (1 = yes).\cr
#' \code{period} \tab final year of case accrual (surrogate for the changing formulations of oral contraceptives over time). \cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Berlin JA, Longnecker MP, Greenland S. Meta-analysis of epidemiologic dose-response data. 
#' Epidemiology. 1993 May 1:218-28. 
#' @keywords data
NULL


#' Six published studies on the relation between alcohol intake and cardiovascular disease risk
#'
#' @name alcohol_cvd
#' @description The dataset reports the summarized dose-response results from six observational
#' studies on the relation between alcohol intake and vascular disease risk. Four are case-control studies,
#' two prospective (cumulative-incidence data).
#' @docType data
#' @format A data frame with 25 observations on the following 8 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author of the study.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted "relative risks".\cr
#' \code{se} \tab standard error for the logarithm of the adjusted "relative risks".\cr
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


#' Four case-control studies on the relation between Body Mass Index and renal cell cancer
#'
#' @name bmi_rc
#' @description The dataset reports the summarized dose-response results from four cases-control studies on the relation
#'  Body Mass Index and renal cell cancer 
#' 
#' @docType data
#' @format A data frame with 33 observations on the following 13 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author and year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{interval} \tab intervals for the categories of bmi.\cr
#' \code{bmi} \tab assigned bmi levels.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{control} \tab number of controls for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{or} \tab adjusted odds ratios for each exposure level.\cr
#' \code{lb_or} \tab lower bound for the confidence limits of the adjusted odds ratios.\cr
#' \code{ub_or} \tab upper bound for the confidence limits of the adjusted odds ratios.\cr
#' \code{logor} \tab natural logarithm of the adjusted odds ratios.\cr
#' \code{se_logor} \tab standard error for the logarithm of the adjusted odds ratios.\cr
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


#' Eight published studies on the relation between alcohol intake and colorectal cancer
#'
#' @name alcohol_crc
#' @description The dataset reports the summarized dose-response results from eight prospective
#' studies on the relation between alcohol intake and colorectal cancer risk.
#' @docType data
#' @format A data frame with 48 observations on the following 7 variables:
#' \tabular{ll}{
#' \code{id} \tab label for author's names (id variable).\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{peryears} \tab amount of person-time for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted "relative risks".\cr
#' \code{se} \tab standard error for the logarithm of the adjusted "relative risks".\cr
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


#' Four published studies on the relation between alcohol intake and lung cancer
#'
#' @name alcohol_lc
#' @description The dataset reports the summarized dose-response results from four prospective
#' studies on the relation between alcohol intake and lunger cancer.
#' @docType data
#' @format A data frame with 20? observations on the following 7 variables:
#' \tabular{ll}{
#' \code{id} \tab label for author's names (id variable).\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{peryears} \tab amount of person-time for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted "relative risks".\cr
#' \code{se} \tab standard error for the logarithm of the adjusted "relative risks".\cr
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
#' 
#' @references
#' Crippa, A., Orsini, N. Dose-response meta-analysis of differences in means. BMC medical research methodology. 2016 Aug 2;16(1):91.
#'
#' @keywords data
NULL


#' Eleven prospective studies on the relation between coffee consumption and stroke risk
#'
#' @name coffee_stroke
#' @description The dataset reports the summarized dose-response results from eleven prospective studies on the relation between coffee consumption
#' and risk of stroke.
#'
#' @docType data
#' @format A data frame with 68 observations on the following 12 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author of the studies.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{rr} \tab adjusted risk estimates for each exposure level.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{nordic} \tab indicator variable for the study to be conducted in the nordic countries (1 = yes).\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Larsson, S. C., Orsini, N. (2011). Coffee consumption and risk of stroke: a dose-response 
#' meta-analysis of prospective studies. American journal of epidemiology, 174(9), 993-1001.
#' 
#' @keywords data
NULL


#' Twenty-one prospective studies on the relation between coffee consumption and all-cause mortality
#'
#' @name coffee_mort
#' @description The dataset reports the summarized dose-response results from twenty-one prospective studies on the relation between coffee consumption
#' and all-cause mortality.
#'
#' @docType data
#' @format A data frame with 109 observations on the following 11 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{gender} \tab factor variable for the gender of the partecipants.\cr
#' \code{area} \tab factor variable for the study location.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Crippa A, Discacciati A, Larsson SC, Wolk A, Orsini N. Coffee Consumption and Mortality from All Causes, Cardiovascular Disease, and Cancer: 
#' A Dose-Response Meta-Analysis. Am J Epidemiol. 2014 Aug 24. pii: kwu194.
#' 
#' @keywords data
NULL


#' Additional two prospective studies on the relation between coffee consumption and all-cause mortality
#'
#' @name coffee_mort_add
#' @description The dataset reports the summarized dose-response results from two additional prospective 
#' studies on the relation between coffee consumption and all-cause mortality. The studies do not report
#' information on the number of cases and participants/person-time.
#'
#' @docType data
#' @format A data frame with 109 observations on the following 11 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{gender} \tab factor variable for the gender of the partecipants.\cr
#' \code{area} \tab factor variable for the study location.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Crippa A, Discacciati A, Larsson SC, Wolk A, Orsini N. Coffee Consumption and Mortality from All Causes, Cardiovascular Disease, and Cancer: 
#' A Dose-Response Meta-Analysis. Am J Epidemiol. 2014 Aug 24. pii: kwu194.
#' 
#' @keywords data
NULL


#' Thirteen prospective studies on the relation between coffee consumption and cardiovascular mortality
#'
#' @name coffee_cvd
#' @description The dataset reports the summarized dose-response results from thirteen prospective studies on the relation between coffee consumption
#' and cardiovascular mortality.
#'
#' @docType data
#' @format A data frame with 100 observations on the following 12 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{gender} \tab character variable for the gender of the partecipants.\cr
#' \code{area} \tab character variable for the study location.\cr
#' \code{smoking} \tab character variable for the type of smoking adjustment.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Crippa A, Discacciati A, Larsson SC, Wolk A, Orsini N. Coffee Consumption and Mortality from All Causes, Cardiovascular Disease, and Cancer: 
#' A Dose-Response Meta-Analysis. Am J Epidemiol. 2014 Aug 24. pii: kwu194.
#' 
#' @keywords data
NULL


#' Eight prospective studies on the relation between coffee consumption and cancer mortality
#'
#' @name coffee_cancer
#' @description The dataset reports the summarized dose-response results from eight prospective studies on the relation between coffee consumption
#' and cancer mortality.
#'
#' @docType data
#' @format A data frame with 59 observations on the following 11 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{gender} \tab factor variable for the gender of the partecipants.\cr
#' \code{area} \tab factor variable for the study location.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Crippa A, Discacciati A, Larsson SC, Wolk A, Orsini N. Coffee Consumption and Mortality from All Causes, Cardiovascular Disease, and Cancer: 
#' A Dose-Response Meta-Analysis. Am J Epidemiol. 2014 Aug 24. pii: kwu194.
#' 
#' @keywords data
NULL


#' Eleven prospective studies on the relation between milk consumption and all-cause mortality
#'
#' @name milk_mort
#' @description The dataset reports the summarized dose-response results from eleven prospective studies on the relation between milk consumption
#' and all-cause mortality.
#'
#' @docType data
#' @format A data frame with 50 observations on the following 12 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{rr} \tab adjusted risk estimates.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Larsson SC, Crippa A, Orsini N, Wolk A, Michaelsson K. Milk consumption and mortality from all causes, 
#' cardiovascular disease, and cancer: a systematic review and meta-analysis. Nutrients. 2015 Sep 11;7(9):7749-63.
#'  
#' @keywords data
NULL


#' Fourteen case-control studies on the relation between alcohol consumption and esophageal cancer
#' @name alcohol_esoph
#' @description The dataset reports the summarized dose-response results from fourteen case-control studies on the relation between alcohol consumption
#' and esophageal squamous cell carcinoma.
#'
#' @docType data
#' @format A data frame with 63 observations on the following 8 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{type} \tab code for study design.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{logrr} \tab natural logarithm of the adjusted odds ratio.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted odds ratio\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Rota M, Bellocco R, Scotti L, Tramacere I, Jenab M, Corrao G, La Vecchia C, Boffetta P, 
#' Bagnardi V. Random-effects meta-regression models for studying nonlinear dose-response 
#' relationship, with an application to alcohol and esophageal squamous cell carcinoma. 
#' Statistics in medicine. 2010 Nov 20;29(26):2679-87.
#'  
#' @keywords data
NULL


#' Six studies on the relation between fish consumption and rheumatoid arthritis risk
#' @name fish_ra
#' @description The dataset reports the summarized dose-response results from six studies on the relation between fish consumption
#' and rheumatoid arthritis risk
#'
#' @docType data
#' @format A data frame with 22 observations on the following 12 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci") or person-years (type = "ir") for each exposure level.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{rr} \tab adjusted risk estimates.\cr
#' \code{lrr} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{urr} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted odds ratio.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted odds ratio\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Di Giuseppe D, Crippa A, Orsini N, Wolk A. Fish consumption and risk of rheumatoid 
#' arthritis: a dose-response meta-analysis. Arthritis research & therapy. 
#' 2014 Sep 30;16(5):446.
#'  
#' @keywords data
NULL


#' Twelve studies on the relation between red meat and bladder cancer
#'
#' @name red_bc
#' @description The dataset reports the summarized dose-response results from twelve studies on the relation between red meat consumption
#' and bladder cancer.
#'
#' @docType data
#' @format A data frame with 74 observations on the following 15 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose0} \tab original assigned dose levels, with unit of measurement defined in the "unit" column.\cr
#' \code{dose} \tab assigned dose levels (converted (if needed) in gm/day).\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci" or "cc") or person-years (type = "ir") for each exposure level.\cr
#' \code{rr} \tab adjusted risk estimates.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{area} \tab geographical area of the published study.\cr
#' \code{unit} \tab unit of measurement for red meat consumption (for dose0).\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Crippa A, Larsson SC, Discacciati A, Wolk A, Orsini N. Red and processed meat consumption and 
#' risk of bladder cancer: a dose-response meta-analysis of epidemiological studies. 
#' European journal of nutrition. 2016 Dec 22:1-3.
#'  
#' @keywords data
NULL


#' Ten studies on the relation between processed meat and bladder cancer
#'
#' @name process_bc
#' @description The dataset reports the summarized dose-response results from ten studies on the relation between processed meat consumption
#' and bladder cancer.
#'
#' @docType data
#' @format A data frame with 73 observations on the following 15 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author.\cr
#' \code{year} \tab year of publication.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose0} \tab original assigned dose levels, with unit of measurement defined in the "unit" column.\cr
#' \code{dose} \tab assigned dose levels (converted (if needed) in gm/day).\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects (type = "ci" or "cc") or person-years (type = "ir") for each exposure level.\cr
#' \code{rr} \tab adjusted risk estimates.\cr
#' \code{lb} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{ub} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' \code{area} \tab geographical area of the published study.\cr
#' \code{unit} \tab unit of measurement for red meat consumption (for dose0).\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Crippa A, Larsson SC, Discacciati A, Wolk A, Orsini N. Red and processed meat consumption and 
#' risk of bladder cancer: a dose-response meta-analysis of epidemiological studies. 
#' European journal of nutrition. 2016 Dec 22:1-3.
#'  
#' @keywords data
NULL


#' Simulated data for one-stage dose-response meta-analysis
#'
#' @name sim_os
#' @description The dataset contains simulated data from 9 case-control studies.
#'
#' @docType data
#' @format A data frame with 27 observations on the following 11 variables:
#' \tabular{ll}{
#' \code{xcati} \tab category limits for the continuous exposure.\cr
#' \code{id} \tab id of the studies.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose levels.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{rr} \tab adjusted risk estimates for each exposure level.\cr
#' \code{lrr} \tab lower bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{urr} \tab upper bound for the confidence limits of the adjusted risk estimates.\cr
#' \code{logrr} \tab natural logarithm of the adjusted risk estimates.\cr
#' \code{se} \tab standard error for the logarithm of the adjusted risk estimates.\cr
#' }
#'
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references
#' 
#' Larsson, S. C., Orsini, N. (2011). Coffee consumption and risk of stroke: a dose-response 
#' meta-analysis of prospective studies. American journal of epidemiology, 174(9), 993-1001.
#' 
#' @keywords data
NULL
