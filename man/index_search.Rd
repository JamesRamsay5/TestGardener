\name{index_search}
\alias{index_search}
\title{
Ensure that estimated score index is global
}
\description{
Multiple minima are found quite often in the data fitting function that is 
minimized using function \code{indexfun}, and in roughly 10 percent of 
the estimates there is a minimum that it lower than that detected.  The
function searches a mesh of 101 points for minima, computes the fitting
function at the minima, and assigns the location of the global minimum as
the replacement index if the location differs by more than 0.5 from 
the value identified by \code{index_fun}.  The function values and their
first two derivatives are also replaced.}
\usage{
  index_search(SfdList, chcemat, index, Fval, DFval, D2Fval, indexind=1:N)
}
\arguments{
  \item{SfdList}{A list vector containing specifications of surprisal
  curves for each item.}
  \item{chcemat}{An \code{N} by \code{n} matrix containing indices of chosen 
  items for each test taker.}
  \item{index}{A vector containing all the score index values.}
  \item{Fval}{A vector containing the \code{N} function values.}
  \item{DFval}{A vector containing the \code{N} first derivative values.}
  \item{D2Fval}{A vector containing the \code{N} second derivative values.}
  \item{indexind}{A vector containing indices of values to be processed.}
}
\value{
A named list object containing objects produced from analyzing the simulations,
one set for each simulation:
  \item{index:}{A vector containing all the score index values
    including those that are altered.}
  \item{Fval:}{A vector containing the N function values
    included those that are altered.}
  \item{DFval:}{A vector containing the N first derivative values
    included those that are altered.}
  \item{D2Fval:}{A vector containing the N second derivative values
    included those that are altered.}
  \item{changeindex:}{Indices of the index values that are altered}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{index_fun}}
}
\examples{
#  Search for values of index that are not at the global minimum of the
#  fitting function and replace them as well as their function and
#  derivative values associated with the fine grid value nearest the 
#  the global minimum.
chcemat <- Quant_13B_problem_chcemat
key     <- Quant_13B_problem_key
SfdList <- Quant_13B_problem_parmList$SfdList
index   <- Quant_13B_problem_parmList$index
Fval    <- Quant_13B_problem_parmList$Fval
DFval   <- Quant_13B_problem_parmList$DFval
D2Fval  <- Quant_13B_problem_parmList$D2Fval
Result  <- index_search(SfdList, chcemat, index, Fval, DFval, D2Fval)
changeindex <- Result$changeindex
print(paste("Number changed =",length(changeindex)))
change  <- index[changeindex] - Result$index[changeindex]
}
