\name{ICC.plot}
\alias{ICC.plot}
\title{Plot probability and surprisal curves for a selection of test or scale 
       items.}
\usage{
  ICC.plot(scrfine, WfdList, dataList, Qvec, 
           binctr = NULL, data_point = FALSE, ci = FALSE, 
           plotType="S", Wrng=c(0,5), DWrng=c(-0.2, 0.2), plotindex=1:n, 
           titlestr = NULL, scopevec = rep(0, length(plotindex)),  
           plotTitle = TRUE,autoplot = FALSE, plotMissing = TRUE, 
           plotrange=c(min(scrfine),max(scrfine)), shaderange = NULL,  
           ttlsz = NULL, axisttl = NULL, axistxt = NULL, 
           lgdlab = NULL, lgdpos = "bottom")
}
\description{
ICC.plot plots each item in argument \code{plotindex} in turn after 
function\code{Wbinsmth()} has used spline smoothing to estimate item and 
option characteristic curves.
}
\arguments{
  \item{scrfine}{A vector of 101 plotting points.}
  \item{WfdList}{A numbered list object produced by a TestGardener analysis of 
    a test.  Its length is equal to the number of items in the test or questions 
    in the scale.  Each member of \code{WfdList} is a named list containing 
    information computed during the analysis.
    These named lists contain these objects:
    \describe{
      \item{scrfine:}{A vector of length 101 containing plotting points.}
      \item{Wfd:}{A functional data object containing the \code{M} surprisal 
        curves for a question.}
      \item{M:}{The number of options.}
      \item{Pbin:}{A matrix containing proportions at each bin.}
      \item{Wbin:}{A matrix containing surprisal values at each bin.}
      \item{indfine:}{A vector of length 101 with equally spaced score values used
      	for plotting.}
      \item{Pmatfine:}{A matrix of probabilities over a fine mesh.}
      \item{Wmatfine:}{A matrix of surprisal values over a fine mesh.}
      \item{DWmatfine:}{A matrix of the values of the first derivative of 
        surprisal curves over fine mesh.}
      \item{D2Wmatfine:}{A matrix of the values of the second derivative of 
        surprisal curves over fine mesh.}
      \item{PStdErr:}{A matrix of the values of the standard error of 
        probability curves at each bin.}
      \item{WStdErr:}{A matrix of the values of the standard error of 
        surprisal curves at each bin.}
    }
  }
  \item{dataList}{A list that contains the objects needed to analyse the test or 
    rating scale with the following fields:
    \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
             N is the number of examinees or respondents and
             n is the number of items.
             Entries in the matrices are the indices of the options chosen.
             Column i of U is expected to contain only the integers 
             \code{1,...,noption}.}  
      \item{optList:}{A list vector containing the numerical score values 
        assigned to the options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is 
       	a numeric vector of length n containing the indices the right answers.  
       	Otherwise, it is NULL.}
      \item{WfdPar:}{An fdPar object for the defining the surprisal curves.}
      \item{noption:}{A numeric vector of length n containing the numbers of 
        options for each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{scrrng:}{A vector of length 2 containing the limits of observed 
        sum scores.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{scrvec:}{A vector of length N containing the examinee or 
        respondent sum scores.}
      \item{itemvec:}{A vector of length n containing the question or item 
        sum scores.}
      \item{percntrnk:}{A vector length N containing the sum score 
        percentile ranks.}
      \item{thetaQnt:}{A numeric vector of length 2*nbin + 1 containing 
        the bin boundaries alternating with the bin centers. These are initially 
        defined as \code{seq(0,100,len=2*nbin+1)}.}
      \item{Wdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 
        5, 25, 50, 75 and 95.}
    }
  }
  \item{Qvec}{A vector of five marker percentile values.  For
    plotting over information, this is replaced by Qvec_al returned
    as parList$Qvec_al.}
  \item{binctr}{A vector of bin center values.  If the plot is over
    arc length or information, \code{binctr} is modified before calling
  \code{Wbinsth.plot} by the command
  \code{binctr_al = pracma::interp1(indfine, alfine, binctr)}, and
    argument \code{binctr} is replaced by \code{binctr_al}.}
  \item{data_point}{A logical value indicating whether to plot the data points.}
  \item{ci}{A logical value indicating whether to plot the confidence limits.}
  \item{plotType}{Type(s) of plot, default as "P" for probability, 
    can also be "W" for surprisal, 
    "DW" for sensitivity, and any combination of the three}
  \item{Wrng}{A vector of length 2 specifing the plotting range for surprisal 
    values.}
  \item{DWrng}{A vector of length 2 specifing the plotting range for 
    sensitivity values.}
  \item{plotindex}{A vector of indices of items to be plotted.}
  \item{titlestr}{plot title}
  \item{scopevec}{A numeric vector of length n containing item scope values.}
  \item{plotTitle}{indicator of showing the plot title, default as TRUE}
  \item{autoplot}{indicator for plotting all items in a batch}
  \item{plotMissing}{Determine if plot the extra option for missing/spoiled responses.}
  \item{plotrange}{A vector of length 2 containing the plot boundaries of
  	 the score index interval.}
  \item{shaderange}{a list of length 2 vector(s); set if users want to gray 
    out specific score range(s)}      
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
  \item{lgdpos}{legend position, could be set as "None" to remove the legend.}
}
\value{
  A list vector is returned which is of the length of argument \code{plotindex}.  
  Each member of the vector is a \code{gg} or \code{ggplot} object for the 
  associated \code{plotindex value}.  Each plot can be displayed using the 
  \code{print} command.
  The plots of item power are produced as a side value even if no output object 
  is specified in the call to the function.
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
  \code{\link{Sensitivity.plot}},
  \code{\link{Power.plot}},
  \code{\link{Entropy.plot}},
  \code{\link{Wbinsmth}}
}
\examples{
#  Example 1.  Display the item surprisal curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
WfdList  <- Quant_13B_problem_parList$WfdList
Qvec     <- Quant_13B_problem_parList$Qvec
binctr   <- Quant_13B_problem_parList$binctr
titlestr <- NULL
#  plot the curves for the first question
scrfine <- seq(0,100,len=101)
oldpar <- par(no.readonly=TRUE)
indfine <- seq(0,100,len=101)
ICC.plot(indfine, WfdList, Quant_13B_problem_dataList, Qvec, binctr, 
         plotType = c("S", "P"), plotindex=6)
par(oldpar)
#  Example 2.  Display the item probability and surprisal curves for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}


