\name{index_fun}
\alias{index_fun}
\title{Compute optimal scores}
\description{
The percentile score index values are estimated for each person.  The estimates
minimize the negative log likelihoods, which are a type of surprisal.  The
main optimization method is a safe-guarded Newton-Raphson method.  

For any iteration the method uses only those scores that are within the interior  
of the interval [0,100] or at a boundary with a first derivative that would
take a step into the interior, and have second derivative values exceeding the 
value of argument \code{crit}.  Consequently the number of values being 
optimized decrease on each iteration, and iterations cease when either 
all values meet the convergence criterion or are optimized on a 
boundary, or when the number of iterations reaches \code{itermax}.
At that point, if there are any interior scores still associated with 
either non-positive second derivatives or values that exceed 
\code{crit}, the minimizing value along a fine mesh is used.

If \code{itdisp} is positive, the number of values to be estimated
are printed for each iteration.
}
\usage{
  index_fun(index, SfdList, chcemat, itermax = 20, crit = 0.001, 
              itdisp = FALSE)
}
\arguments{
  \item{index}{A vector of size \code{N} containing initial values for score 
    indices in the interval [0,100].}
  \item{SfdList}{A list vector of length equal to the number of questions.  
    Each member contains eight results for the surprisal curves 
    associated with a question.}
  \item{chcemat}{A matrix number of rows equal to the number of examinees or 
    respondents, and number of columns equal to number of items.  The 
    values in the matrix are indices of choices made by each respondent
    to each question.}
  \item{itermax}{Maximum number of iterations for computing the optimal
    index values.  Default is 20.}
  \item{crit}{Criterion for convergence of optimization.  Default is 1e-8.}
  \item{itdisp}{If TRchcematE, results are displayed for each iteration.}
}
\value{
A named list with these members:
  \item{index_out:}{A vector of optimized score index value.} 
  \item{Fval:}{The negative log likelihood criterion.}
  \item{DFval:}{The first derivative of the negative likelihood.} 
  \item{D2Fval:}{The second derivative of the negative likelihood.} 
  \item{iter:}{The number iterations used.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{index_distn}},
  \code{\link{Ffun}},
  \code{\link{DFfun}},
  \code{\link{index2info}},
  \code{\link{scoreDensity}}
}
\examples{
  #  Optimize the indices defining the data fits for the first five examinees
  #  input the choice indices in the 1000 by 24 choice index matrix
  chcemat   <- Quant_13B_problem_chcemat
  #  First set up the list object for surprisal curves computed from 
  #  initial index estimates.
  SfdList   <- Quant_13B_problem_dataList$SfdList
  #  Their initial values are the percent rank values ranging over [0,100]
  index_in  <- Quant_13B_problem_dataList$percntrnk[1:5]
  #  set up choice indices for first five examinees
  chcemat_in <- chcemat[1:5,]
  #  optimize the initial indices
  indexfunList <- index_fun(index_in, SfdList, chcemat_in)
  #  optimal index values
  index_out    <- indexfunList$index_out
  #  The surprisal data fit values
  Fval_out     <- indexfunList$Fval
  #  The surprisal data fit first derivative values
  DFval_out    <- indexfunList$DFval
  #  The surprisal data fit second derivative values
  D2Fval_out   <- indexfunList$D2Fval
  #  The number of index values that have not reached the convergence criterion
  active_out   <- indexfunList$active
}
