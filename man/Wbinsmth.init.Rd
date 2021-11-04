\name{Wbinsmth.init}
\alias{Wbinsmth.init}
\title{
Set up initial coefficient values for computing surprisal curves.
}
\description{
Excellent approximations of the $K$ by $M-1$ coefficient matrix for an item
can be computed using ordinary spline smoothing, where $K$ is the number
of basis functions used to compute surprisal curves and $M$ is the number
of options for a particular item.  This function is called
by function \code{make.dataList.R} prior to the analysis of the data.
}
\usage{
Wbinsmth.init(percntrnk, nbin, WfdPar, grbg, optList, U)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{percntrnk}{Sum scores are first converted to ranks and then the ranks are 
  normalized to be within the interval 0 to 100.}
  \item{nbin}{The number of bins used to compute bin proportions.}
  \item{WfdPar}{A functional parameter object defining a set of surprisal curves
  and optionally defining a roughness penalty.}
  \item{grbg}{A vector of indices for each item of the category used to 
  record missing or illegal response.  If there are no such responses, 
  this will be NULL.}
  \item{optList}{A list vector containing for each item an item label, 
  labels for each option for that item, and the scores assigned to 
  each option.  Labels are optional, but option scores must be
  provided.}
  \item{U}{A matrix of choice index values for each test taker and each item.}
}
\value{
An list vector \code{WfdList} containing for each item a functional data
object named \code{Wfd} that is defined by the initial coefficient 
approximations.  The approximate coefficient matrix for item 
\code{i} is in \code{WfdList[[i]]$Wfd$coefs}.}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.

\url{http://testgardener.azurewebsites.net}
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Wbinsmth}},
  \code{\link{Wbinsmth.plot}},
  \code{\link{thetafun}}
}