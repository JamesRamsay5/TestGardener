\name{Sensitivity.plot}
\alias{Sensitivity.plot}
\title{
  Plots all the sensitivity curves for selected items or questions.}
\description{
  A sensitivity curve for an option is the first derivative of the corresponding
  surprisal curve.  Its values can be positive or negative, and the size of the departure
  from zero at any point on the curve is the amount information contributed by that curve 
  to locating the value of an examinee or respondent on the score index continuum.
}
\usage{
Sensitivity.plot(scrfine, WfdList, Qvec, dataList, plotindex=1:n, 
                 plotrange=c(min(scrfine),max(scrfine)),   
                 key=NULL, titlestr=NULL, saveplot=FALSE, width=c(-0.2,0.2), 
                 ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
}
\arguments{
  \item{scrfine}{A vector of length nfine (usually 101) containing equally spaced
  	points spanning the plotrange.  Used for plotting.}
  \item{WfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  Its length
    is equal to the number of items in the test or questions in the scale.  
    Each member of \code{WfdList} is a named list containing information computed during the analysis.
    These named lists contain these objects:
    \describe{
      \item{Wfd:}{A functional data object containing the \code{M} surprisal curves
        for a question.}
      \item{M:}{The number of options.}
      \item{Pbin:}{A matrix containing proportions at each bin.}
      \item{Wbin:}{A matrix containing surprisal values at each bin.}
      \item{Pmatfine:}{A matrix of probabilities over a fine mesh.}
      \item{Wmatfine:}{A matrix of surprisal values over a fine mesh.}
      \item{DWmatfine:}{A matrix of the values of the first derivative of 
        surprisal curves over fine mesh.}
      \item{D2Wmatfine:}{A matrix of the values of the second derivative of 
        surprisal curves over fine mesh.}
    }
  }
  \item{Qvec}{The values of the five marker percentiles.}
  \item{dataList}{A list that contains the objects needed to analyse the test or rating scale
  	with the following fields:
    \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
        N is the number of examinees or respondents and n is the number of items.
        Entries in the matrices are the indices of the options chosen.
        Column i of U is expected to contain only the integers 1,...,noption.}  
      \item{optList:}{A list vector containing the numerical score values assigned to the
        options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is a numeric vector 
       	of length n containing the indices the right answers.  Otherwise, it is NULL.}
      \item{WfdPar:}{An fdPar object for the defining the surprisal curves.}
      \item{noption:}{A numeric vector of length n containing the numbers of options for 
      	each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{scrrng:}{A vector of length 2 containing the limits of observed sum scores.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{scrvec:}{A vector of length N containing the examinee or respondent sum scores.}
      \item{itemvec:}{A vector of length n containing the question or item sum scores.}
      \item{percntrnk:}{A vector length N containing the sum score percentile ranks.}
      \item{thetaQnt:}{A numeric vector of length 2*nbin + 1 containing the bin boundaries
      	alternating with the bin centers. These are initially defined as 
      	seq(0,100,len=2*nbin+1).}
      \item{Wdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 5, 25, 50, 75 and 95.}
    }
  }
  \item{plotindex}{A set of integers specifying the numbers of the items or questions 
    to be displayed.}
  \item{plotrange}{A vector of length 2 containing the plot boundaries within or
    over the score index interval c(0,100).}
  \item{key}{A integer vector of indices of right answers.  
    If the data are rating scales, this can be NULL.}
  \item{titlestr}{A title string for plots.}
  \item{saveplot}{A logical value indicating whether the plot should be saved to a pdf file.}
  \item{width}{A vector of length 2 defining the lower and upper limits on the ordinate 
    for the plots.}
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
}
\details{
  Sensitivity curves for each question indexed in the \code{index} argument.  
  A request for a keystroke is made for each question.
}
\value{
  A list vector is returned which is of the length of argument \code{plotindex}.  
  Each member of the vector is a \code{gg} or \code{ggplot} object for the associated
  \code{plotindex value}.  Each plot can be displayed using the \code{print} command.
  The plots of item power are produced as a side value even if no output object is
  specified in the call to the function.
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
  \code{\link{Power.plot}},
  \code{\link{Entropy.plot}},
  \code{\link{ICC.plot}}
}
\examples{
#  Example 1.  Display the option sensitivity curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
WfdList <- Quantshort_parList$WfdList
Qvec    <- Quantshort_parList$Qvec
scrfine <- seq(0,100,len=101)
Sensitivity.plot(scrfine, WfdList, Qvec, Quantshort_dataList, 
                 plotindex=1)
#  Example 2.  Display the option sensitivity curves for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
}