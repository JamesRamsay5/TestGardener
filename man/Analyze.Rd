\name{Analyze}
\alias{Analyze}
\title{
  Analyze test or rating scale data defined in \code{dataList}.
}
\description{
   The test or rating scale data have already been processed by function \code{make_dataList} or
   or other code to produce the list object \code{dataList}.  The user defines a list vector
   \code{ParameterList} which stores results from a set of cycles of estimating surprisal curves
   followed by estimating optimal score index values for each examinee or respondent.  These
   score index values are within the interval [0,100].  The number of analysis cycles is the
   length of the \code{parmList} list vector.
}
\usage{
  Analyze(index, indexQnt, dataList, ncycle=10, itdisp=FALSE, verbose=FALSE) 
}
\arguments{
  \item{index}{A vector of \code{N} score index values for the examinees or respondents.
  These values are in the percent interval [0,100].}
  \item{indexQnt}{A vector of length \code{2*nbin + 1} where \code{nbin} is the 
  number of bins containing score index values.  The vector begins with the lower
  boundary 0 and ends with the upper boundary 100.  In between it alternates 
  between the bin center value and the boundary separating the next bin.}
  \item{dataList}{A list that contains the objects needed to analyse the test or rating scale
  	with the following fields:
    \describe{
      \item{chce:}{A matrix of response data with N rows and n columns where
             N is the number of examinees or respondents and
             n is the number of items.
             Entries in the matrices are the indices of the options chosen.
             Column i of chce is expected to contain only the integers 1,...,noption.}  
      \item{optList:}{A list vector containing the numerical score values assigned to the
       options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is a numeric vector 
       	of length n containing the indices the right answers.  Otherwise, it is NULL.}
      \item{Sfd:}{An fd object for the defining the surprisal curves.}
      \item{noption:}{A numeric vector of length n containing the numbers of options for 
      	each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{scrrng:}{A vector of length 2 containing the limits of observed sum scores.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{scrvec:}{A vector of length N containing the examinee or respondent sum scores.}
      \item{itemvec:}{A vector of length n containing the question or item sum scores.}
      \item{percntrnk:}{A vector length N containing the sum score percentile ranks.}
      \item{indexQnt:}{A numeric vector of length 2*nbin + 1 containing the bin boundaries
      	alternating with the bin centers. These are initially defined as 
      	seq(0,100,len=2*nbin+1).}
      \item{Sdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 5, 25, 50, 75 and 95.}
    }
  }
  \item{ncycle}{The number of cycles executed by \code{Analyze}.}
  \item{itdisp}{If TRUE, the progress of the iterations within each cycle for estimating
  	index are reported.}
  \item{verbose}{If TRUE, the stages of analysis within each cycle for estimating
  	index are reported.}
}
\details{
   The cycling process is described in detail in the references, and displayed in R code 
   in the vignette \code{SweSATQuantitativeAnalysis}.
}
\value{
  The list vector \code{parmList} where each member is a named list object containing 
  the results of an analysis cycle.  These results are:
  \item{index:}{The optimal estimates of the score index values for the 
    examinees/respondents.  This is a vector of length \code{N}.}
  \item{indexQnt:}{A vector of length 2*nbin+1 containing bin boundaries alternating
    with bin edges.}
  \item{SfdList:}{A list vector containing results from the estimation of surprisal 
    curves. The list vector is of length \code{n}, the number of questions or items in 
    the test of rating scale.  For details concerning these results, see function 
    \code{Sbinsmth()}.}
  \item{meanF:}{For each person, the mean of the optimal fitting function values.}
  \item{binctr:}{A vector of length \code{nbin} containing the bin centers within 
    the interval [0,100].}
  \item{bdry:}{A vector of length \code{nbin+1} containing the bin boundaries.}
  \item{freq:}{A vector of length \code{nbin} containing the number of score 
    index values in the bins.  An score index value is within a bin if it is less than 
    or equal to the upper boundary and greater than the lower boundary.  The first 
    boundary also contains zero values.}
  \item{pdf_fd}{Functional probability curves}
  \item{logdensfd:}{A functional data object defining the estimate of the log of the 
    probability density function for the distribution of the score index values.}
  \item{C:}{The normalizing value for probability density functions.  
    A density value is computed by dividing the exponential of the log density value 
    by this constant.}
  \item{denscdf:}{The values over a fine mesh of the cumulative probability 
    distribution function. These values start at 0 and end with 1 and are increasing.  
    Ties are often found at the upper boundary, so that using these values for 
    interpolation purposes may require using the vector \code{unique(denscdf)}.}
  \item{indcdf}{Equally spaced index values to match the number in denscdf.}
  \item{Qvec}{Locations of the marker percents.}
  \item{index}{The positions of each test taker on the score index continuum.}
  \item{Fval:}{A vector of length N containing the values of the negative log 
    likelihood fitting criterion.}
  \item{DFval:}{A vector of length N containing the values of the first
    derivative of the negative log likelihood fitting criterion.}
  \item{D2Fval:}{A vector of length N containing the values of the second
    derivative of the negative log likelihood fitting criterion.}
  \item{active:}{A vector of length N of the activity status of the 
    values of index.  If convergence was not achieved, the value is TRUE,
    otherwise FALSE.}  
  \item{infoSurp:}{The length of the space curve defined by the surprisal curves.}
}
\references{
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
  Journal of Educational and Behavioral Statistics, 45, 297-315.
  
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
  information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{make_dataList},
        \link{TG_analysis},
        \link{index_distn},
        \link{index2info},
        \link{index_fun},
        \link{Sbinsmth}}
}

\examples{
\dontrun{
  #  Example 1:  Input choice data and key for the short version of the 
  #  SweSAT quantitative multiple choice test with 24 items and 1000 examinees
  #  input the choice data as 1000 strings of length 24
  #  setup the input data list object
  dataList <- Quant_13B_problem_dataList
  #  define the initial examinee indices and bin locations
  index    <- dataList$percntrnk
  indexQnt <- dataList$indexQnt
  #  Set the number of cycles (default 10 but here 5)
  ncycle <- 5
  parmListvec <- Analyze(index, indexQnt, ncycle=ncycle, dataList, 
                         verbose=TRUE) 
  #  two column matrix containing the mean fit and arclength values
  #  for each cycle
  HALsave <- matrix(0,ncycle,2)
  for (icycle in 1:ncycle) {
    HALsave[icycle,1] <- parmListvec[[icycle]]$meanF
    HALsave[icycle,2] <- parmListvec[[icycle]]$infoSurp
  }
  #  plot the progress over the cycles of mean fit and arc length
  par(mfrow=c(2,1))
  plot(1:ncycle, HALsave[,1], type="b", lwd=2, 
       xlab="Cycle Number",ylab="Mean H")
  plot(1:ncycle, HALsave[,2], type="b", lwd=2, 
       xlab="Cycle Number", ylab="Arc Length")}
}
