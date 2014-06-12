\name{generate.simulation}
\alias{generate.simulation}
\title{Generate Simulation Data}
\description{This function generates a simulated data set using a specific setting.}
\usage{
generate.simulation(pc1 = c(1,0.5), pc2 = c(-0.5,1), weight = c(1,0.5),
  rho=0.4, rho1=0.2, err.var=2, beta=c(0.5,-0.5,0.5,-0.3,0.6,-0.3),
  n=1000)
}
\arguments{
  \item{pc1}{first eigen-direction (or principal component direction) of covariates}
  \item{pc2}{second eigen-direction (or principal component direction) of covariates}
  \item{weight}{the "true" single index direction}
  \item{rho}{correlation within \code{U} or \code{V}}
  \item{rho1}{correlation between \code{U} and \code{V}}
  \item{err.var}{noise of quadratic response surface}
  \item{beta}{"true" coefficients of quadratic response surface}
  \item{n}{sample size}
  }
\value{
  \item{data}{the simulated data set}
  \item{cov_mat}{covariance matrix of \code{U} and \code{V}}
  \item{rho_mat}{correlation matrix of \code{U} and \code{V}}
}
\author{Vivian W. Chan, Huan Cheng, John L. Michela, Mu Zhu}
\examples{
\dontrun{
junk=generate.simulation()
junk$cov_mat
junk$rho_mat
pairs(junk$data)
}
}
