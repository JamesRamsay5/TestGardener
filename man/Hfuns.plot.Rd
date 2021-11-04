\name{Hfuns.plot}
\alias{Hfuns.plot}
\title{
Plot a selection of fit criterion H functions and their first two derivatives.
}
\description{
These plots indicate whether an appropriate minimum of the fitting criterion was found.
The value of theta should be at the function minimum, the first derivative be close to
zero there, and the second derivative should be positive.  If these conditions are not
met, it may be worthwhile to use function \code{thetafun} initialized with an approximate  
minimum value of score index theta to re-estimate the value of theta.}
\usage{
  Hfuns.plot(theta, WfdList, U, plotindex=1)
}
\arguments{
  \item{theta}{The entire vector of estimated values of theta.}
  \item{WfdList}{The list vector of length \code{n} containing the 
  estimated surprisal curves.}
  \item{U}{The entire N by n matrix of choice indices.}
  \item{plotindex}{A subset of the integers \code{1:N}.}
}
\details{
The curves are displayed in three vertically organized panels along with 
values of theta and the values and first two derivative values of the 
fit criterion.  If more than one index value is used, a press of the Enter 
or Return key moves to the next index value.
}
\value{
A list vector is returned which is of the length of argument \code{plotindex}.  Each
member of the vector is a \code{gg} or \code{ggplot} object for the associated
\code{plotindex value}.  Each plot can be displayed using the \code{print} command.
The plots of item power are produced as a side value even if no output object is
specified in the call to the function.
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
  \code{\link{thetafun}},
  \code{\link{Hfun}},
  \code{\link{DHfun}}
}
\examples{
#  Example 1.  Display fit criterion values and derivatives for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
index <- 1
theta    <- Quant_parList$theta
WfdList  <- Quant_parList$WfdList
U        <- Quant_dataList$U
Hval     <-  Hfun(theta[index], WfdList, U[index,])
DHResult <- DHfun(theta[index], WfdList, U[index,])
DHval    <- DHResult$DH
D2Hval   <- DHResult$D2H
Hfuns.plot(theta, WfdList, U, plotindex=index)
#  Example 2.  Compute score index values theta for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
}
