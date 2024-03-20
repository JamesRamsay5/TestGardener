\name{Scope_plot}
\alias{Scope_plot}
\title{
Plot the score index \code{index} as a function of arc length.
}
\description{Arc length or scope is the distance along the space curved traced 
  out as score index \code{index} increases from 0 to 100.  
  It is measured in bits and is remains 
  unchanged if the score index continuum is modified.}
\usage{
  Scope_plot(infoSurp, infoSurpvec, titlestr=NULL)
}
\arguments{
  \item{infoSurp}{This is a positive real number indicating the total length of the 
  space curve.  It is expressed in terms of numbers of bits.}
  \item{infoSurpvec}{A vector of length 101 containing equally-spaced arc-length distances  
  along the test information curve.}
  \item{titlestr}{A string for the title of the data.}
}
\value{
  A \code{gg} or \code{ggplot} object defining the plot of infoSurp along the test information 
  curve as a function of the score index \code{index}.  This is displayed by the \code{print} command. 
  The plot is automatically displayed as a side value even if no return object is specified in 
  the calling statement.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{index2info}}
}
\examples{
#
#  Example 1.  Display the arc length curve for the 
#  SweSAT multiple choice test with 24 items and 1000  examinees
#
infoSurpvec <- Quant_13B_problem_infoList$infoSurpvec
infoSurp    <- Quant_13B_problem_infoList$infoSurp
oldpar <- par(no.readonly=TRUE)
Scope_plot(infoSurp, infoSurpvec)
on.exit(oldpar)
}

