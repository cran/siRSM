\name{surface.stats}
\alias{surface.stats}
\alias{surface.stats.main}
\title{Statistics of Quadratic Response Surface}
\description{Computes stationary point, 1st and 2nd principal axes of the quadratic response surface, shapes and curvatures along the congruence and incongruence lines, as well as various parameters useful for \code{plot.siRSM}. Mostly for internal use only.}
\usage{
surface.stats(obj)
surface.stats.main(b, xlim, ylim)
}
\arguments{
  \item{obj}{an object of class \code{siRSM}, often the result of \code{siRSM}}
  \item{b}{coefficients of the response surface, e.g., \code{obj$coef} where \code{obj} is an object of class \code{siRSM}}
  \item{xlim}{x-boundaries in 3D (x,y,z)-perspective plot --- used when \code{surface.stats.main} is called internally by \code{plot.siRSM}}
  \item{ylim}{y-boundaries in 3D (x,y,z)-perspective plot --- used when \code{surface.stats.main} is called internally by \code{plot.siRSM}}
}
\value{
 \item{u0}{x-coordinates of surface's stationary point}
 \item{v0}{y-coordinates of surface's stationary point}
 \item{p10}{intercept in xy-plane of surface's 1st principal axis}
 \item{p11}{slope in xy-plane of surface's 1st principal axis}
 \item{p20}{intercept in xy-plane of surface's 2nd principal axis}
 \item{p21}{slope in xy-plane of surface's 2nd principal axis}
 \item{ax.congr}{slope of surface along the congruence line, \code{u-v=0}}
 \item{ax2.congr}{curvature of surface along the congruence line, \code{u-v=0}}
 \item{ax.incongr}{slope of surface along the incongruence line, \code{u+v=0}}
 \item{ax2.incongr}{curvature of surface along the incongruence line, \code{u+v=0}}
 \item{pl}{intersection of 1st axis at lower border, to be used by \code{plot.siRSM}}
 \item{ph}{intersection of 1st axis at upper border, to be used by \code{plot.siRSM}}
 \item{sl}{intersection of 2nd axis at lower border, to be used by \code{plot.siRSM}}
 \item{sh}{intersection of 2nd axis at upper border, to be used by \code{plot.siRSM}}
}
\references{
Edwards JR (2002), "Alternatives to difference scores: Polynomial regression analysis and response surface methodology," in \cite{Measuring and Analyzing Behavior in Organizations: Advances in Measurement and Data Analysis} (F. Drasgow, N. Schmitt, Eds.), Jossey-Bass, pp. 350-400.
}
\author{Huan Cheng, Mu Zhu}
