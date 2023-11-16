\name{ICC}
\alias{ICC}
\title{
  Objects required for representing a single test or scale item
}
\description{
  This generates an S3 object that contains for a single item a surprisal fd
  object, the probability and surprisal binned data values, probability and 
  surprisal values of a mesh of 101 score index values and probability and 
  surprisal standard errors for each bin.  The objects in the argument list 
  are computed from option choices by \code{smooth.ICC.}}
\usage{
  ICC(Item, M, Sfd, Pbin, Sbin, Pmatfine, Sarrayfine, PStdErr, SStdErr, 
      ItemArcLen, itemStr, optStr)
}
\arguments{
  \item{Item}{An integer indicating the position of the item within
        1,...,n}
  \item{Sfd}{A functional data object containing the \code{M} 
        surprisal curves for a single item.}
  \item{M}{The number of options in a single item.}
  \item{Pbin}{A matrix containing proportions at each bin.}
  \item{Sbin}{A matrix containing surprisal values at each bin.}
  \item{Pmatfine}{A matrix of probabilities over a fine mesh.}
  \item{Sarrayfine}{An 101 by M by 3 array  of surprisal values over 
      a fine mesh of 101 values in \code{0,100}. The layers contain the 
      surprisal values, and their first and second derivatives, respectively.}
  \item{PStdErr}{An nbin by M  matrix of the values of the standard error 
      of probability estimates.}
  \item{SStdErr}{An nbin by M  matrix of the values of the standard error
      of surprisal estimates.}
  \item{ItemArcLen}{A non-negative surprisal value for the arc length of
      the item information curve.}
  \item{itemStr}{A title string for the item.}
  \item{optStr}{A character vector containing title strings for options with 
  each item.}
}
\details{This object replaces one of the corresponding objects stored in  
list object \code{SfdList}.}
\value{An S3 object with class ICC.}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{ICC_plot}}
}
\examples{
#  Example 1.  Display the item probability and surprisal curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  estimate the initial surprisal curves
index     <- Quant_13B_problem_dataList$percntrnk
indexQnt  <- Quant_13B_problem_parmList$indexQnt
Qvec      <- Quant_13B_problem_parmList$Qvec
SfdResult <- Sbinsmth(index, Quant_13B_problem_dataList)
# set up item 1
SfdList1  <- SfdResult$SfdList[[1]]
Item      <- SfdList1$Item
M         <- SfdList1$M
Sfd       <- SfdList1$Sfd
n         <- length(Quant_13B_problem_dataList)
# set up item 1
}

