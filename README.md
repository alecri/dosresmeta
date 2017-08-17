## Multivariate Dose-Response Meta-Analysis (dosresmeta)

The package `dosresmeta` consists of a collection of functions to estimate dose-response relations from 
summarized dose-response data for both continuous and binary outcomes, and to combine them according to 
principles of (multivariate) random-effects model. The methodology is illustrated in the referenced article.


### Info on the `dosresmeta` package

The package is available on the Comprehensive R Archive Network (CRAN), with info at the related web page (https://CRAN.R-project.org/package=dosresmeta). 
A development website is available on GitHub (https://github.com/alecri/dosresmeta).

For a short summary of the package, refer to the main help page by typing:

```r
help("dosresmeta-package")
```

in R after installation (see below). 

### Installation

The last version officially released on CRAN can be installed directly within R by typing:

```r
install.packages("dosresmeta")
```

A version still under developement is avaiable on [GitHub](https://github.com/alecri/dosresmeta) and can be installed by typing:

```r
install.packages("devtools")
devtools::install_github("alecri/dosresmeta")
```

### R code in published articles

Several peer-reviewed articles and documents provide R code illustrating methodological developments or replicating 
substantive results. 
An updated version of the code can be found at the GitHub (https://github.com/alecri) or personal web page 
(https://alecri.github.io/software/dosresmeta.html) of the package maintainer.

### References:

Crippa A, Orsini N. Multivariate Dose-Response Meta-Analysis: the dosresmeta R Package. 
*Journal of Statistical Software*, Code Snippets,. 2016; 72(1), 1-15. doi:10.18637/jss.v072.c01. [freely available [here](https://alecri.github.io/downloads/jss1256.pdf)]

Crippa A, Orsini N. Dose-response meta-analysis of differences in means. *BMC Medical Research Methodology*. 2016 Aug 2;16(1):91. [freely available [here](https://www.researchgate.net/publication/305804878_Dose-response_meta-analysis_of_differences_in_means)] [GitHub repository at this [link](https://github.com/alecri/differences-in-mean)]

Discacciati A, Crippa A, Orsini N. Goodness of fit tools for dose-response meta-analysis of binary outcomes. *Research Synthesis Methods*. 2015 Jan 1. doi: 10.1002/jrsm.1194. [freely available [here](http://onlinelibrary.wiley.com/doi/10.1002/jrsm.1194/pdf)] [GitHub repository at this [link](https://github.com/anddis/goodness-of-fit-meta-analysis)]




