\name{ci.index}
\alias{ci.index}
\title{Condidence Intervals for Single Index}
\description{Computes 95\% bootstrap confidence intervals for the single index.}
\usage{
ci.index(y, U, V, B=100, use.parallel=TRUE, ...)
}
\arguments{
  \item{y}{vector, response}
  \item{U}{matrix, whose columns are covariates for factor one}
  \item{V}{matrix, whose columns are covariates for factor two}
  \item{B}{number of bootstap samples to take}
  \item{use.parallel}{if TRUE, exploits multiple cores by using \pkg{foreach}, \pkg{doSNOW}, etc}
  \item{...}{other arguments for \code{siRSM}, e.g., \code{interaction.only=TRUE}, \code{trial=10}}
}
\value{A data.frame specifying the lower 2.5\% and upper 97.5\% confidence limits, the mean, and the standard error for each coordinate of the index.
}
\author{Mu Zhu}
\keyword{bootstrap}
\keyword{parallel computing}
\keyword{foreach}