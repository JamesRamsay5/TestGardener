\name{make.dataList}
\alias{make.dataList}
\title{
  Make a list object containing information required for analysis of choice data.
}
\description{
  The returned named list object contains 11 named members, which are described 
  in the value section below.
}
\usage{
make.dataList(U, key, optList, grbg=rep(0,n), scrrng=NULL, titlestr=NULL,  
              nbin=nbinDefault(N),NumBasis=NULL,  WfdPar=NULL,
              jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
              quadwrd=FALSE, verbose=FALSE)
}
\arguments{
  \item{U}{A matrix with rows corresponding to examinees or respondents, 
    and columns to questions or items.}
  \item{key}{If the data are a multiple choice test with only weights 0 and 1, 
    a vector of length n containing the indices of the right answers.  
    Otherwise, NULL.}
  \item{optList}{A list vector of length number of questions. A member 
    contains a vector of score values assigned to each answer or option by the 
    test designer.}
  \item{grbg}{A vector of length indicating which option for each item contains 
    missing or illegal choice values.  If 0, there is no such option.}
  \item{scrrng}{A vector of length two containing the initial and final values  
    for the interval over which test scores are to be plotted.  Default is 
    minimum and maximum sum score.}
  \item{titlestr}{A string to be used as a title in plots and other displays.}
  \item{nbin}{The number of bins for containing proportions of examinees 
    choosing options.  The default is computed by a function that uses the 
    number of examinees.}
  \item{NumBasis}{The number of spline basis functions used to represent 
    surprisal curves. The default is computed by a function that uses the 
    number of examinees.}
  \item{WfdPar}{A functional parameter object specifying a basis, a linear 
    differential operator, a smoothing parameter, a boolian constant
    indicating estimation, and a penalty matrix.  The default object defines
    an order 5 spline with operator \code{Lfd = 3} and smoothing parameter
    \code{1e4}. The number of basis functions depends on the size of N.}
  \item{jitterwrd}{A boolian constant: TRUE implies adding a small random value 
    to each sum score value prior to computing percent rank values.}
  \item{PcntMarkers}{Used in plots of curves to display marker or reference 
    percentage points for abscissa values in plots.}
  \item{quadwrd}{If TRUE a default order 5 basis is constructed.
    If FALSE the order and number of of basis functions are 3 so as to
    support the use of the second derivative in function \code{thetafun}.}
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
  are involved.  A minimum of 7 and maximum of 24 basis functions is assigned 
  by default.  The number of basis function can be set by the user or otherwise
  by a simple algorithm  used to adjust the number of bins to the number 
  \code{N} of examinees.
  
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
        N is the number of examinees or respondents and n is the number of items.
        Entries in the matrices are the indices of the options chosen.
        Column i of U is expected to contain only the integers \code{1,...,noption}.}  
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
N <- dim(Quantshort_U)[1]
n <- dim(Quantshort_U)[2]

optList <- list()
grbg = matrix(0,n,1)
for (item in 1:n){
  noptioni <- length(unique(Quantshort_U[,item]))
  scorei <- rep(0,noptioni)
  scorei[Quantshort_key[item]] <- 1
  optList[[item]] <- scorei
  grbg[item] <- length(scorei)
}
optList <- list(itemLab=NULL, optLab=NULL, optScr=optList)

#  Set up the dataList object containing the objects necessary
#  for further display and analyses

Quantshort_dataList <- make.dataList(Quantshort_U, Quantshort_key, optList, grbg)

#  Example 2:  Input choice data and key for the Symptom Distress Scale 
#  with 13 items and 473 examinees.
#  input the choice data as 473 strings of length 13
N <- dim(SDS_U)[1]
n <- dim(SDS_U)[2]
# --------- Define the option score values for each item ---------
optList <- list()
grbg <- matrix(0,n,1)
for (item in 1:n){
  scorei <- c(0:4,0)
  optList[[item]] <- scorei
  grbg[item] <- length(scorei)
}
optList <- list(itemLab=NULL, optLab=NULL, optScr=optList)
#  largest observed sum score is 37 
scrrng <- c(0,37)
SDS_dataList <- make.dataList(SDS_U, SDS_key, optList, grbg, scrrng)
}
