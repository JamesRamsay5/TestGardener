\name{ConditionalSimulation}
\alias{ConditionalSimulation}
\title{
  Simulation-based Estimation of Efficiency of a TestGardener Model
}
\description{
  Given a psychometric model defined by a WfdCell object, a set of simulated samples are constructed and then analyzed using the model.  The efficiency of the estimates of each of 101 equally spaced values of the score index theta over [0,100] is assessed by computing the root-mean-squared errors and biasses at each of these values.  A graphical display of the results is generated.
}
\usage{
ConditionalSimulation(dataList, parList, theta.pop=seq(0,100,len=101), nsample=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataList}{A list object containing required information about the data structure.
  }
  \item{parList}{A list object containing required information about the model structure.}
  \item{theta.pop}{A set of score index values at which efficiency is to be assessed. 
  Defaults to the 101 values 0, 1, ..., 100.}
  \item{nsample}{The number of simulated samples.  Defaults to 1000.}
}
\value{
  A named list with named members as follows:
  \describe{
    \item{sumscrsave:}{A 101 by nsample matrix of sum scores.}
    \item{thetasave:}{A 101 by nsample matrix of estimated score index values.}
    \item{musave:}{A 101 by nsample matrix of estimated expected sum score values.}
    \item{alsave:}{A 101 by nsample matrix of estimated arc length or information values.}
    \item{theta.pop:}{The true score index values.}
    \item{mu.pop:}{The true expected sum score values.}
    \item{al.pop:}{The true arc length or information values.}
    \item{Qvec:}{The model marker percentages.}
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
  \code{\link{scorePerformance}},
  \code{\link{Usimulate}}
}
