\name{siRSM}
\alias{siRSM}
\alias{siRSM.default}
\alias{plot.siRSM}
\alias{plot.siRSM.default}
\alias{print.siRSM}
\alias{print.siRSM.default}
\title{Single-index Response Surface Model}
\description{This is the main function for users to call.}
\usage{
\method{siRSM}{default}(y, U, V, trial, interaction.only=FALSE, use.parallel=TRUE)

\method{print}{siRSM}(x, ...)

\method{plot}{siRSM}(x, ...)
}
\arguments{
  \item{y}{response vector}
  \item{U}{matrix, whose columns are covariates for first factor}
  \item{V}{matrix, whose columns are covariates for second factor, must have same dimension as \code{U}}
  \item{trial}{number of different initial values to try --- if missing, defaults to \code{2*K}, where \code{K = ncol(U) = ncol(V)}}
  \item{interaction.only}{fit an interaction-only (rather than a full-quadratic) model}
  \item{use.parallel}{if TRUE, tells underlying utility function \code{multi.run} to exploit multiple cores by using \pkg{foreach}, \pkg{doSNOW}, etc}
  \item{x}{an object of class \code{siRSM}, often result of \code{siRSM}}
  \item{...}{arguments for \code{plot} include:
  \code{xname}=character string, name for first (composite) factor;
  \code{yname}=character string, name for second (composite) factor;
  \code{zname}=character string, name for response;
  \code{center}=if 'zero', centers the surface plot at \code{(0,0)}}
}
\value{An object of class \code{siRSM}. For most users, the most useful elements are:
 \item{w}{the estimated single index vector}
 \item{coef}{the coefficients of the corresponding response surface}
}
\author{Huan Cheng, Mu Zhu}
\keyword{single index models}
\keyword{response surface models}
\keyword{parallel computing}