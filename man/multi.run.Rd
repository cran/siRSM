\name{multi.run}
\alias{multi.run}
\title{Fit siRSM by Trying Multiple Initial Values}
\description{Internal function called by \code{siRSM}. Do NOT use.}
\usage{
multi.run(y, X, Z, rep, interaction.only=FALSE, use.parallel=TRUE)
}
\arguments{
  \item{y}{vector, response}
  \item{X}{matrix, whose columns are covariates for factor one}
  \item{Z}{matrix, whose columns are covariates for factor two, must be of same size as X}
  \item{rep}{number of different initial values to try --- if unspecified, the default is twice the dimension of the index}
  \item{interaction.only}{fit an interaction-only model}
  \item{use.parallel}{if TRUE, exploits multiple cores by using \pkg{foreach}, \pkg{doSNOW}, etc}
}
\note{In these more rudimentary functions (not to be called by user), \code{U} is referred to as \code{X}, and \code{V} as \code{Z}, due to "historical reasons". Currently, multiple RANDOM initial values are used. In the future, we'd like to switch to using a space-filling design.}
\author{Huan Cheng, Mu Zhu}
\keyword{parallel computing}
\keyword{foreach}