\name{siRSM-package}
\alias{siRSM-package}
\docType{package}
\title{siRSM: Single-index Response Surface Models}
\description{This package provides tools for fitting what we call "single-index response surface models", that is, models of the form \code{y=f(u,v)+e}, where \code{f(u,v) ~ 1 + u + v + u^2 + uv + v^2}, \code{u=t(w)\%*\%U}, and \code{v=t(w)\%*\%V}.}
\details{\tabular{ll}{
Package: \tab siRSM\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-07-15\cr
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
# load illustrative data set provided within the package
data(EdwardsMBA)

# parse the variables
y=EdwardsMBA[,'AVGSAT4']
U=EdwardsMBA[,c('PCPRE','DMPRE','EIPRE','MRPRE')]
V=EdwardsMBA[,c('PCACT','DMACT','EIACT','MRACT')]

# fit the model
m1=siRSM(y,U,V)

# look at the model
m1
plot(m1)

# F-test of surface curvature
surface.test(m1)

# inference on the index (using just 10 bootstrap samples here)
ci.index(y,U,V,B=10)

# obtain statistics of the response surface (conditional on the estimated index)
surface.stats(m1)

# inference for these surface statistics (again, conditional on the estimated index)
ci.surface(m1)

# fit an interaction-only model
m2=siRSM(y,U,V,interaction.only=TRUE)
plot(m2)
}
}
