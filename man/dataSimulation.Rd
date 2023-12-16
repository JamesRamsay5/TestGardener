\name{dataSimulation}
\alias{dataSimulation}
\title{Simulation Based Estimates of Error Variation of Score Index Estimates}
\description{Estimate sum score,s score index values index and test information 
values bias and mean squared errors using simulated data.}
\usage{
  dataSimulation(dataList, parmList, nsample = 1000)
}
\arguments{
  \item{dataList}{The list object set up by function \code{make_dataList}.}
  \item{parmList}{The list object containing objects computed by function
  \code{Analyze}.}
  \item{nsample}{The number of simulated samples.}
}
\value{
A named list object containing objects produced from analyzing the simulations,
one set for each simulation:
\item{sumscr:}{Sum score estimates}
\item{index:}{Score index estimates}
\item{mu:}{Expected sum score estimates}
\item{info:}{Total arc length estimates}
\item{index.pop:}{True or population score index values}
\item{mu.pop:}{Expected sum score population values}
\item{info.pop:}{Total test length population values}
\item{n:}{Number of items}
\item{nindex:}{Number of index values}
\item{indfine:}{Fine mesh over score index range} 
\item{Qvec:}{Five marker percentages: 5, 25, 50, 75 and 95}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{scorePerformance}}
}
