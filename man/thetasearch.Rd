\name{thetasearch}
\alias{thetasearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ensure that estimated score index is global
}
\description{
Multiple minima are found quite often in the data fitting function that is 
minimized using function \code{thetafun}, and in roughly 10 percent of 
the estimates there is a minimum that it lower than that detected.  The
function searches a mesh of 101 points for minima, computes the fitting
function at the minima, and assigns the location of the global minimum as
the replacement theta if the location differs by more than 0.5 from 
the value identified by \code{thetafun}.  The function values and their
first two derivatives are also replaced.}
\usage{
  thetasearch(WfdList, U, theta, Hval, DHval, D2Hval, thetaind=1:N)
}
\arguments{
  \item{WfdList}{A list vector containing specifications of surprisal
  curves for each item.}
  \item{U}{An N by n matrix containing indices of chosen items for each
  test taker.}
  \item{theta}{A vector containing all the score index values.}
  \item{Hval}{A vector containing the N function values.}
  \item{DHval}{A vector containing the N first derivative values.}
  \item{D2Hval}{A vector containing the N second derivative values.}
  \item{thetaind}{A vector containing indices of values to be processed.}
}
\value{
A named list object containing objects produced from analyzing the simulations,
one set for each simulation:
  \describe{
    \item{theta:}{A vector containing all the score index values
    including those that are altered.}
    \item{Hval:}{A vector containing the N function values
    included those that are altered.}
    \item{DHval:}{A vector containing the N first derivative values
    included those that are altered.}
    \item{D2Hval:}{A vector containing the N second derivative values
    included those that are altered.}
    \item{changeindex:}{Indices of the theta values that are altered}
  }
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{thetafun}}
}
\examples{
#  Search for values of theta that are not at the global minimum of the
#  fitting function and replace them as well as their function and
#  derivative values associated with the fine grid value nearest the 
#  the global minimum.
WfdList <- Quant_13B_problem_parList$WfdList
theta   <- Quant_13B_problem_parList$theta
Hval    <- Quant_13B_problem_parList$Hval
DHval   <- Quant_13B_problem_parList$DHval
D2Hval  <- Quant_13B_problem_parList$D2Hval
Result  <- thetasearch(WfdList, Quant_13B_problem_U, theta, 
                       Hval, DHval, D2Hval, thetaind=1:1000)
changeindex <- Result$changeindex
print(paste("Number changed =",length(changeindex)))
change  <- theta[changeindex] - Result$theta[changeindex]
oldpar <- par(no.readonly=TRUE)
hist(change,101)
par(oldpar)
}
