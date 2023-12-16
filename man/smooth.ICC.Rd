\name{smooth.ICC}
\alias{smooth.ICC}
\title{
  Smooth binned probability and surprisal values to make an \code{ICC} object.
}
\description{
An N by n matrix of positive integer choice index values is transformed to
an nbin by M matrix of probability values by iteravely minimizing the sum 
of squared errors for bin values.}
\usage{
smooth.ICC(x, item, index, dataList, indexQnt=seq(0,100, len=2*nbin+1), 
                       wtvec=matrix(1,n,1), iterlim=20, conv=1e-4, dbglev=0)
}
\arguments{
  \item{x}{An ICC object}
  \item{item}{Index of item being set up.}
  \item{index}{A vector of length N containing score index values for each
  person.}
  \item{dataList}{A list object set up by function \code{make.dataList}
  containing objects set up prior to an analysis of the data.}
  \item{indexQnt}{A vector of length 2*nbin + 1 containing, in sequence,
  the lower boundary of a bin, its midpoint, and the upper boundary.}
  \item{wtvec}{A vector of length n containing wseights for items.}
  \item{iterlim}{An integer specifying the maximum number of optimizations.}
  \item{conv}{A convergence criterion a little larger than 0.}
  \item{dbglev}{One of integers 0 (no optimization information),
  1 (one line per optimization) or 2 (complete optimization display).}
}
\value{An S3 class ICC object for a single item.}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\examples{
# example code to be set up
}