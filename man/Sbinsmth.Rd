\name{Sbinsmth}
\alias{Sbinsmth}
\title{Estimate the option probability and surprisal curves.}
\description{The surprisal curves for each item are fit to the surprisal 
transforms of choice probabilities for each of a set of bins of current
performance values \code{index}.  The error sums of squares are minimized
by the surprisal optimization \code{smooth.surp} in the \code{fda}
package.  The output is a list vector of length \code{n} containing the
functional data objects defining the curves.}
\usage{
  Sbinsmth(index, dataList, indexQnt=seq(0,100, len=2*nbin+1), 
           wtvec=matrix(1,n,1), iterlim=20, conv=1e-4, dbglev=0)
}
\arguments{
  \item{index}{A vector of length N containing current values of score index 
    percentile values.}  
  \item{dataList}{A list that contains the objects needed to analyse the test or 
    rating scale.}
  \item{indexQnt}{A vector of length 2*n+1 containing the sequence of bin boundary and 
    bin centre values.}
  \item{wtvec}{A vector of length \code{n} of weights on observations.  Defaults to
    all ones.}
  \item{iterlim}{The maximum number of iterations used in optimizing surprisal
    curves.  Defaults to 20.}
  \item{conv}{Convergence tolerance.  Defaults to 0.0001.}
  \item{dbglev}{Level of output within \code{Sbinsmth}.  If 0, no output,
    if 1 the error sum of squares and slope on each iterations, and if 2
    or higher, results for each line search iteration with function \code{lnsrch}.}
}
\details{
  The function first bins the data in order to achieve rapid estimation of the 
  option surprisal curves. The argument \code{indexQnt} contains the sequence 
  of bin boundaries separated by the bin centers, so that it is of length 
  \code{2*nbin + 1} where \code{nbin}
  is the number of bins.  
  These bin values are distributed over the percentile interval
  [0,100] so that the lowest boundary is 0 and highest 100.  
  Prior to the call to \code{Sbinsmth} these boundaries are computed so that 
  the numbers of values of \code{index} falling in the bins are roughly equal.  
  It is important that the number of bins be chosen so that the bins contain 
  at least about 25 values.  
  
  After the values of \code{index} are binned, the proportions that the bins 
  are chosen for each question and each option are computed.  Proportions of 
  zero are given NA values.
  
  The positive proportions are then converted to surprisal values where 
  surprisal = -log_M (proportion) where log_M is the logarithm with base M, 
  the number of options associated with a question.  
  Bins with zero proportions are assigned a surprisal that is appropriately 
  large in the sense of being in the range of the larger surprisal values 
  associated with small but positive proportions. This surprisal value
  is usually about 4.
  
  The next step is to fit the surprisal values for each question by a 
  functional data object that is smooth, passes as closely as possible to an
  option's surprisal values, and has values consistent with being a surprisal 
  value.  The function \code{smooth.surp()} is used for this purpose.  The
  arc length of thme item information curve is also computed.
  
  Finally the curves and other results for each question are saved in object
  \code{SfdList}, a list vector of length n, and the list vector is returned.
}
\value{
  The optimized numbered list object \code{SfdList} with length \code{n} 
  that provides data on the probability and surprisal data and curves.  
  The 12 objects for each item are as follows:
  \item{Sfd:}{A surprisal functional data object that is used for plotting.
  It also contains the coefficient matrix and functional data basis 
  that define the object.}
  \item{M:}{The number of options, including if needed a final
  option which is for the missing and illegitimate responses.}
  \item{Pbin:}{A \code{nbin} by \code{M} matrix of proportions of 
  choice for each option.}
  \item{Sbin:}{A \code{nbin} by \code{M} matrix of surprisal values
  for each option..}
  \item{indfine:}{A fine mesh of 101 equally spaced score index values
  over the interval [0,1].}
  \item{Pmatfine:}{A 101 by \code{M} matrix of probability values
  at each of the fine mesh points \code{indfine}.}
  \item{Smatfine:}{A 101 by \code{M} matrix of surprisal values
  at each of the fine mesh points \code{indfine}.}
  \item{DSmatfine:}{A 101 by \code{M} matrix of surprisal 
  first derivative values at each of the fine mesh points 
  \code{indfine}.}
  \item{D2Smatfine:}{A 101 by \code{M} matrix of surprisal 
  second derivative values at each of the fine mesh points 
  \code{indfine}.}
  \item{PSrsErr:}{The standard error for probability over
  the fine mesh.}
  \item{PSrsErr:}{The standard error for surprisal over
  the fine mesh.}
  \item{itemScope:}{The length of the item info curve.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{ICC_plot}},
  \code{\link{Sbinsmth}}
}
\examples{
#  Example 1.  Display the initial probability and surprisal curves for the 
#  first item in the short SweSAT multiple choice test with 24 items and 
#  1000 examinees.
#  Note: The scope is 0 at this point because it is computed later 
#  in the analysis.
dataList <- Quant_13B_problem_dataList
index    <- dataList$percntrnk
#  Carry out the surprisal smoothing operation
SfdResult   <- Sbinsmth(index, dataList)
  \dontrun{
  #  Set up the list object for the estimated surprisal curves
  SfdList     <- SfdResult$SfdList
  #  The five marker percentage locations for (5, 25, 50, 75, 95)
  binctr      <- dataList$binctr
  Qvec        <- dataList$PcntMarkers
  #  plot the curves for the first question
  scrfine   <- seq(0,100,len=101)
  ICC_plot(scrfine, SfdList, dataList, Qvec, binctr,
            data_point = TRUE, plotType = c("S", "P"), 
            Srng=c(0,3), plotindex=1)
  }
}
