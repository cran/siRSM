\name{ci.surface}
\alias{ci.surface}
\title{Conditional Inference of Response Surface}
\description{Computes 95\% bootstrap confidence intervals for various features of the quadratic response surface, CONDITIONAL on the single index.}
\usage{
ci.surface(obj, B=500, use.parallel=TRUE)
}
\arguments{
  \item{obj}{an object of class \code{siRSM}, typically result from \code{siRSM}}
  \item{B}{number of bootstrap samples}
  \item{use.parallel}{if TRUE, exploits multiple cores by using \pkg{foreach}, \pkg{doSNOW}, etc}
}
\value{In all components below, the lower 2.5\% and upper 97.5\% confidence limits, the mean, and the standard error are given:
  \item{stationary.point}{stationary point, (\code{u0, v0})}
  \item{prin.ax.1}{1st principal axis, intercept (\code{p10}) and slope (\code{p11})}
  \item{prin.ax.2}{2nd principal axis, intercept (\code{p20}) and slope (\code{p21})}
  \item{beta}{coefficients defining the quadratic response surface, \code{b0, b1, ..., b5}}
  \item{line.congr}{slope (\code{ax}) and curvature (\code{ax2}) along the congruence line, \code{u-v=0}}
  \item{line.incongr}{slope (\code{ax}) and curvature (\code{ax2}) along the congruence line, \code{u+v=0}}
}
\author{Huan Cheng, Mu Zhu}
\keyword{bootstrap}
\keyword{parallel computing}
\keyword{foreach}