\name{siRSM-package}
\alias{siRSM-package}
\docType{package}
\title{siRSM: Single-index Response Surface Models}
\description{This package provides tools for fitting what we call "single-index response surface models", that is, models of the form \code{y=f(u,v)+e}, where \code{f(u,v) ~ 1 + u + v + u^2 + uv + v^2}, \code{u=t(w)\%*\%U}, and \code{v=t(w)\%*\%V}.}
\details{
\tabular{ll}{
Package: \tab siRSM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-12\cr
License: \tab GPL-2\cr
}}
\author{Huan Cheng, Mu Zhu

Maintainer: Mu Zhu <m3zhu@uwaterloo.ca>
}
\references{
Cheng H, Zhu M, Chan VW, Michela, JL (in preparation), "Single-index response surface models".
}
\keyword{single index models}
\keyword{response surface models}
\keyword{parallel computing}
\seealso{
\code{\link{siRSM}}, \code{\link{surface.test}}, \code{\link{surface.stats}}, \code{\link{ci.index}}, \code{\link{ci.surface}}
}
\examples{
\dontrun{
# generate some data (using a function provided within this package)
tmp=generate.simulation()
y=tmp$data[,'intent']
U=tmp$data[,c('s1','s2')]
V=tmp$data[,c('e1','e2')]

# fit a single-index RSM
model=siRSM(y,U,V,use.parallel=TRUE)
model

# test the curvature of the quadratic response surface CONDITIONAL on the index
surface.test(model)

# obtain various feature parameters for the quadratic response surface
surface.stats(model)

# plot the surface
plot(model)

# confidence intervals for the index (using just 10 bootstrap samples here)
indstat=ci.index(y,U,V,B=10,use.parallel=TRUE)
indstat

# (conditional) confidence intervals for various surface statistics
surfstat=ci.surface(model)
surfstat
}
}
