\name{Ffuns_plot}
\alias{Ffuns_plot}
\title{
Plot a selection of fit criterion F functions and their first two derivatives.
}
\description{
These plots indicate whether an appropriate minimum of the fitting criterion was found.
The value of index should be at the function minimum, the first derivative be close to
zero there, and the second derivative should be positive.  If these conditions are not
met, it may be worthwhile to use function \code{indexfun} initialized with an approximate  
minimum value of score index index to re-estimate the value of index.}
\usage{
  Ffuns_plot(evalarg, index, SfdList, chcemat, plotindex=1)
}
\arguments{
  \item{evalarg}{A vector containingg the sore index values to be evaluated.}
  \item{index}{The vector of of length \code{N} of score index values.}
  \item{SfdList}{The list vector of length \code{n} containing the 
    estimated surprisal curves.}
  \item{chcemat}{The entire \code{N} by \code{n} matrix of choice indices.}
  \item{plotindex}{A subset of the integers \code{1:N}.}
}
\details{
  The curves are displayed in three vertically organized panels along with 
  values of index and the values and first two derivative values of the 
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
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{index_fun}},
  \code{\link{Ffun}},
  \code{\link{DFfun}}
}
\examples{
#  Example 1.  Display fit criterion values and derivatives for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
chcemat   <- Quant_13B_problem_dataList$chcemat
index     <- Quant_13B_problem_parmList$index
SfdList   <- Quant_13B_problem_parmList$SfdList
plotindex <- 1:3
indfine   <- seq(0,100,len=101)
Ffuns_plot(indfine, index, SfdList, chcemat, plotindex)
}
