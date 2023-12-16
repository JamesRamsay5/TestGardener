\name{Fcurve}
\alias{Fcurve}
\title{
  Construct grid of 101 values of the fitting function
}
\description{
A fast grid of values of the fitting function or one of its first two 
derivatives is constructed for use in function \code{indexsearch}.}
\usage{
  Fcurve(SfdList, chcevec, nderiv=0)
}
\arguments{
  \item{SfdList}{A list vector containing specifications of surprisal
  curves for each item.}
  \item{chcevec}{A N by n matrix containing indices of chosen items for each
  test taker.}
  \item{nderiv}{Integer 0, 1 or 2 to indicate which level of derivative to use.}
}
\value{
A vector of length 101 containing grid values of a derivative of the fitting
function}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{index_search}}
}
\examples{
#  Compute a grid of values of the objective function for locating each  
#  examinee or respondent for the 24-item short form of the SweSAT   
#  quantitative test on the percentile score index continuum [0,100].
chcemat <- Quant_13B_problem_dataList$chcemat
SfdList <- Quant_13B_problem_parmList$SfdList
index   <- Quant_13B_problem_parmList$index
n       <- ncol(chcemat)
#  Fitting function for the first examinee
j <- 1
chcevec   <- as.numeric(chcemat[j,])
Fcurve1   <- Fcurve(SfdList, chcevec, 0)
#  First derivative of the fitting function for the first examinee
DFcurve1  <- Fcurve(SfdList, chcevec, 1)
#  Second derivative of the fitting function for the first examinee
D2Fcurve1 <- Fcurve(SfdList, chcevec, 2)
oldpar  <- par(no.readonly=TRUE)
par(mfrow=c(3,1))
indfine <- seq(0,100,len=101)
plot(indfine, Fcurve1,   type="l", xlab="", ylab="Fitting curve",
     main="Examinee 1")
plot(indfine, DFcurve1,  type="l", xlab="", ylab="First derivative")
points(index[1], 0, pch="o")
abline(0,0,lty=2)
plot(indfine, D2Fcurve1, type="l", 
     xlab="Score index", ylab="Second derivative")
abline(0,0,lty=2)
points(index[1], 0, pch="o")
par(oldpar)
}
