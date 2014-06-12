\name{single.run}
\alias{single.run}
\title{Fit siRSM Once Based on One Initial Value}
\description{Internal function called by \code{siRSM}. Do NOT use.}
\usage{
single.run(data, w0, int.only=FALSE, eps=0.005, precision=1e-6, max.iter=10000, 
  debug=FALSE, trace=TRUE)
}
\arguments{
  \item{data}{\code{y, X, Z} together}
  \item{w0}{initial guess}
  \item{int.only}{fit an interaction-only model}
  \item{eps}{step size for gradient descent}
  \item{precision}{convergence criterion}
  \item{max.iter}{maximum number of iteractions}
  \item{debug}{can be turned on for debugging}
  \item{trace}{if TRUE, reports progress as function runs}
}
\author{Huan Cheng, Mu Zhu}
\note{
In these more rudimentary functions (not to be called by user), \code{U} is referred to as \code{X}, and \code{V} as \code{Z}, due to "historical reasons".
}
