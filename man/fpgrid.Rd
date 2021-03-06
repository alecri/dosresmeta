% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additional.R
\name{fpgrid}
\alias{fpgrid}
\title{Grid with combinations of p for two-order fractional polynomials}
\usage{
fpgrid(p = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3))
}
\arguments{
\item{p}{a numeric vector with the coefficient to be combined.}
}
\value{
A data.frame with the different combinations of p.
}
\description{
Computes the different combinations of p usefull for evaluating two-order fractional polynomials.
}
\examples{

grd <- fpgrid()
head(grd)

}
\references{
Royston, Patrick, and Douglas G. Altman. "Regression using fractional polynomials of 
continuous covariates: parsimonious parametric modelling." Applied Statistics (1994): 429-467.
}
\author{
Alessio Crippa, \email{alessio.crippa@ki.se}
}
