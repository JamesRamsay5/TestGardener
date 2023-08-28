\name{SimulateData}
\alias{SimulateData}
\title{
Simulate Choice Data from a Previous Analysis
}
\description{
Simulation of data using a previous analysis requires only an ICC vector and two objects computed by function \code{theta.distn} along with a specification of the number of simulated the simulated persons.}
\usage{
SimulateData(nsim, indfine, denscdf, WfdList)
}
\arguments{
  \item{nsim}{Number of persons having simulated choices.}
  \item{indfine}{The score index values within [0,100] that are
    associated with the cumulative probability values in 
  \code{denscdf}.}
  \item{denscdf}{The cumulative probability values within [0,1].
    The values have to be discrete, begin with 0 and end with 1.}
  \item{WfdList}{List vector of length n of list vectors for
        item objects.}
}
\details{
Arguments \code{indfine} and \code{denscdf} can be obtained from
the original analysis, but also can be specified to describe
a different distribution of score index values.}
\value{
An nsim by n matrix of integers including 1 and 2 that specify
each person's option choice for each item.}
\references{
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
  Journal of Educational and Behavioral Statistics, 45, 297-315.
  
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
  information-based psychometrics.  Psych, 2, 347-360.
  
  http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{dataSimulation}},
  \code{\link{Usimulate}}
}
\examples{
# example code to be set up 
}