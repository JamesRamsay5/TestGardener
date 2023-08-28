\name{make.dataList}
\alias{make.dataList}
\title{
  Make a list object containing information required for analysis of choice 
  data.
}
\description{
  The list object \code{dataList} contains 22 objects that supply all of the 
  information required to analyze the data.  
  Initial values of the score indices in object \code{theta} and the bin 
  boundaries and centres in object \code{thetaQnt}.
  The returned named list object contains 22 named members, which are described 
  in the value section below.
}
\usage{
make.dataList(U, key, optList, scrrng=NULL, titlestr=NULL,
              nbin=nbinDefault(N), NumBasis=7, WfdPar=NULL,
              jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
              verbose=FALSE)
}
\arguments{
  \item{U}{A matrix with rows corresponding to examinees or respondents, 
    and columns to questions or items.}
  \item{key}{If the data are a multiple choice test with only weights 0 and 1, 
    a vector of length n containing the indices of the right answers.  
    Otherwise, key = NULL, which is the default value. This is the case
    for rating scale questions.}
  \item{optList}{A list vector of length number of questions. A member 
    contains a numeric vector of score values assigned to each answer or option  
    by the test designer. For multiple choice items with 0-1 scroing, the 
    zero values must be included in each item's vector. Scores are not to
    be included for illegal or missing responses, since these are automatically
    added.}
  \item{scrrng}{A numeric vector of length two containing the initial and   
    final values for the interval over which test scores are to be plotted.   
    Default is minimum and maximum sum score.}
  \item{titlestr}{A string to be used as a title in plots and other displays.}
  \item{nbin}{The number of bins for containing proportions of examinees 
    choosing options.  The default is computed by a function that uses the 
    number of examinees.}
  \item{NumBasis}{The number of spline basis functions used to represent 
    surprisal curves. The default is computed by a function that uses the 
    number of examinees.  
    }
  \item{WfdPar}{A functional parameter object specifying a basis, a linear 
    differential operator, a smoothing parameter, a boolian constant
    indicating estimation, and a penalty matrix.  The default object defines
    an order 5 spline with operator \code{Lfd = 3} and smoothing parameter
    \code{1e4}. The number of basis functions depends on the size of N.}
  \item{jitterwrd}{A boolian constant: TRUE implies adding a small random value 
    to each sum score value prior to computing percent rank values.}
  \item{PcntMarkers}{Used in plots of curves to display marker or reference 
    percentage points for abscissa values in plots.}
  \item{verbose}{If TRUE details of calculations are displayed.}
}
\details{
  The score range defined \code{scrrng} should contain all of the sum score 
  values, but can go beyond them if desired.  For example, it may be that no 
  examinee gets a zero sum score, but for reporting and display purposes using 
  zero as the lower limit seems desirable.
  
  The number of bins is chosen so that a minimum of at least about 25 initial 
  percentage ranks fall within a bin.  For larger samples, the number per bin 
  is also larger, making the proportions of choice more accurate.  The number 
  bins can be set by the user, or by a simple algorithm used to adjust the 
  number of bins to the number \code{N} or examinees.
  
  The number of spline basis functions used to represent a surprisal curve 
  should be small for small sample sizes, but can be larger when larger samples 
  are involved.  
  
  There must be at least two basis functions, corresponding to two 
  straight lines.  The norder of this simple spline would
  not exceed 1, corresponding to taking only a single derivative of 
  the resulting spline.  But this rule is bent here to allow higher
  higher derivatives, which will autmatically have values of zero, in 
  order to allow these simple linear basis functions to be used.  This 
  permits direct comparisons of TestGardener models with the many classic 
  item response models that use two or less parameters per item response 
  curve.
    
  Adding a small value to discrete values before computing ranks is considered 
  a useful way of avoiding any biasses that might arise from the way the data 
  are stored.  The small values used leave the rounded jittered values fixed, 
  but break up ties for sum scores.
  
  It can be helpful to see in a plot where special marker percentages 
  5, 25, 50, 75 and 95 percent of the interval [0,100] 
  are located.  The median abscissa value is at 50 per cent for initial 
  percent rank values, for example, but may not be located at the center of   
  the interval after iterations of the analysis cycle.
}
\value{
  A named list with named members as follows:
  \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
        N is number of examinees or respondents and n is number of items.
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
       \item{grbg:}{A logical vector of length number of questions. 
       TRUE for an item indicates that a garbage option must be added to the  
       score values, and FALSE indicates that there are no illegal or missing 
       responses and the number of options is equal to number of score values.}
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
  \code{\link{Analyze}}
}
\examples{
#  Example 1:  Input choice data and key for the short version of the 
#  SweSAT quantitative multiple choice test with 24 items and 1000 examinees
#  input the choice data as 1000 strings of length 24
N <- dim(Quant_13B_problem_U)[1]
n <- dim(Quant_13B_problem_U)[2]
noption <- rep(0,n)
for (i in 1:n) noption[i]  <- 
    length(unique(Quant_13B_problem_U[,i]))
ScoreList <- list() # option scores
for (item in 1:n){
  scorei <- rep(0,noption[item])
  scorei[Quant_13B_problem_key[item]] <- 1
  ScoreList[[item]] <- scorei
}
optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

#  Set up the dataList object containing the objects necessary
#  for further display and analyses

Quant_13B_problem_dataList <- 
    make.dataList(Quant_13B_problem_U, Quant_13B_problem_key, optList)

#  Example 2:  Input choice data and key for the Symptom Distress Scale 
#  with 13 items and 473 examinees.
#  input the choice data as 473 strings of length 13
N <- dim(SDS_U)[1]
n <- dim(SDS_U)[2]
# --------- Define the option score values for each item ---------
optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)
#  largest observed sum score is 37 
scrrng <- c(0,37)
SDS_dataList <- make.dataList(SDS_U, SDS_key, optList, scrrng)
}
