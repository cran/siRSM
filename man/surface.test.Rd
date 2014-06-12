\name{surface.test}
\alias{surface.test}
\title{F-test of Curvature for the Response Surface}
\description{A standard nested F-test of linear (\code{f(u,v)~u+v}), interaction-only (\code{f(u,v)~u+v+I(u*v)}), and full-quadratic (\code{f(u,v)~u+v+I(u^2)+I(u*v)+I(v^2)}) models.}
\usage{
surface.test(object)
}
\arguments{
  \item{object}{an object of class '\code{siRSM}', often result from \code{siRSM}}
}
\author{Huan Cheng, Mu Zhu}
