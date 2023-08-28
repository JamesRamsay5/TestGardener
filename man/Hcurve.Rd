\name{Hcurve}
\alias{Hcurve}
\title{
  Construct grid of 101 values of the fitting function
}
\description{
A fast grid of values of the fitting function or one of its first two 
derivatives is constructed for use in function \code{thetasearch}.}
\usage{
  Hcurve(WfdList, U, nderiv=0)
}
\arguments{
  \item{WfdList}{A list vector containing specifications of surprisal
  curves for each item.}
  \item{U}{A N by n matrix containing indices of chosen items for each
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

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{thetasearch}}
}
\examples{
#  Compute a grid of values of the objective function for locating each  
#  examinee or respondent for the 24-item short form of the SweSAT quantitative  
# test on the percentile score index continuum [0,100].
WfdList <- Quant_13B_problem_parList$WfdList
theta   <- Quant_13B_problem_parList$theta
U       <- Quant_13B_problem_dataList$U
n       <- ncol(U)
#  Fitting function for the first examinee
j <- 1
Umat    <- matrix(U[j,],1,n)
Hcurve1   <- Hcurve(WfdList, Umat, 0)
#  First derivative of the fitting function for the first examinee
DHcurve1  <- Hcurve(WfdList, Umat, 1)
#  Second derivative of the fitting function for the first examinee
D2Hcurve1 <- Hcurve(WfdList, Umat, 2)
oldpar <- par(no.readonly=TRUE)
par(mfrow=c(3,1))
indfine <- seq(0,100,len=101)
plot(indfine, Hcurve1,   type="l", xlab="", ylab="Fitting curve",
     main="Examinee 1")
plot(indfine, DHcurve1,  type="l", xlab="", ylab="First derivative")
points(theta[1], 0, pch="o")
abline(0,0,lty=2)
plot(indfine, D2Hcurve1, type="l", xlab="Score index", ylab="Second derivative")
abline(0,0,lty=2)
points(theta[1], 0, pch="o")
par(oldpar)
}
