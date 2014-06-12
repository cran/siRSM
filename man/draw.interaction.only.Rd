\name{draw.full.quadratic}
\alias{draw.interaction.only}
\alias{draw.full.quadratic}
\title{Functions for plotting various siRSMs}
\description{Internal function called by \code{plot.siRSM}. Do NOT use.}
\usage{
draw.interaction.only(x, xname=NULL, yname=NULL, zname=NULL)
draw.full.quadratic(x, xname=NULL, yname=NULL, zname=NULL, center='scale', 
  debug=FALSE)
}
\arguments{
  \item{x}{an object of class \code{siRSM}}
  \item{xname}{character string, name for first (composite) factor}
  \item{yname}{character string, name for second (composite) factor}
  \item{zname}{character string, name for response}
  \item{center}{if 'scale', centers the surface plot at \code{(0,0)}; otherwise, centers ...}
  \item{debug}{can be turned on for debugging purposes}
}
\author{Huan Cheng, Mu Zhu}
