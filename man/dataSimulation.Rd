\name{dataSimulation}
\alias{dataSimulation}
\title{
Simulation Based Estimates of Root-mean-squared-err of Theta Estimates}
\description{Estimate sum score,s score index values theta and test information 
values bias and mean squared errors using simulated data.}
\usage{
  dataSimulation(dataList, parList, theta.pop = seq(0, 100, len = 101), 
                 nsample = 1000)
}
\arguments{
  \item{dataList}{The list object set up by function \code{make_dataList}.}
  \item{parList}{The list object containing objects compuated by function
  \code{Analyze}.}
  \item{theta.pop}{A vector containing true values of theta to be estimated
  using simulated data.}
  \item{nsample}{The number of simulated samples.}
}
\value{
A named list object containing objects produced from analyzing the simulations,
one set for each simulation:
  \describe{
    \item{sumscr:}{Sum score estimates}
    \item{theta:}{Score index estimates}
    \item{mu:}{Expected sum score estimates}
    \item{al:}{Total arc length estimates}
    \item{thepop:}{True or population score index values}
    \item{mupop:}{Expected sum score population values}
    \item{alpop:}{Total test length population values}
    \item{n:}{Number of items}
    \item{ntheta:}{Number of theta values}
    \item{indfine:}{Fine mesh over score index range} 
    \item{Qvec:}{Five marker percentages: 5, 25, 50, 75 and 95}
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
  \code{\link{scorePerformance}}
}
