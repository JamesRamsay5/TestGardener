\name{thetafun}
\alias{thetafun}
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
  thetafun(theta, WfdList, U, itermax = 20, crit = 0.001, 
           itdisp = FALSE)
}
\arguments{
  \item{theta}{A vector of size \code{N} containing initial values for score 
    indices in the interval [0,100].}
  \item{WfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  
    Its length is equal to the number of items in the test or questions in the 
    scale.  Each member of \code{WfdList} is a named list containing information 
    computed during the analysis.
    The named list in each member contains these objects:
    \itemize{
      \item{Wfd:} {A functional data object containing the \code{M} 
        surprisal curves for a question.}
      \item{M:} {The number of options.}
      \item{Pbin:} {A matrix containing proportions at each bin.}
      \item{Wbin:} {A matrix containing surprisal values at each bin.}
      \item{Pmatfine:} {A matrix of probabilities over a fine mesh.}
      \item{Wmatfine:} {A matrix of surprisal values over a fine mesh.}
      \item{DWmatfine:} {A matrix of the values of the first derivative of 
        surprisal curves over fine mesh.}
      \item{D2Wmatfine:} {A matrix of the values of the second derivative of 
        surprisal curves over fine mesh.}
    }
  }
  \item{U}{A matrix number of rows equal to the number of examinees or 
    respondents, and number of columns equal to number of items.  The 
    values in the matrix are indices of choices made by each respondent
    to each question.}
  \item{itermax}{Maximum number of iterations for computing the optimal
    theta values.  Default is 20.}
  \item{crit}{Criterion for convergence of optimization.  Default is 1e-8.}
  \item{itdisp}{If TRUE, results are displayed for each iteration.}
}
\value{
A named list with these members:
  \describe{
    \item{theta_out:}{A vector of optimized score index value.} 
    \item{Hval:}{The negative log likelihood criterion.}
    \item{DHval:}{The first derivative of the negative likelihood.} 
    \item{D2Hval:}{The second derivative of the negative likelihood.} 
    \item{iter:}{The number iterations used.}
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
  \code{\link{theta.distn}},
  \code{\link{Hfun}},
  \code{\link{DHfun}},
  \code{\link{theta2arclen}},
  \code{\link{testscore}},
  \code{\link{scoreDensity}}
}
