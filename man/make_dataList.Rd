\name{make_dataList}
\alias{make_dataList}
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
  make_dataList(chcemat, scoreList, noption, sumscr_rng=NULL, 
                titlestr=NULL, itemlabvec=NULL, optlabList=NULL,
                nbin=nbinDefault(N), NumBasis=7, jitterwrd=TRUE, 
                PcntMarkers=c( 5, 25, 50, 75, 95), verbose=FALSE)
}
\arguments{
  \item{chcemat}{An N by n matrix.  Column \code{i} must contain the integers 
    from 1 to \code{M_i}, where \code{M_i} is the number of options
    for item \code{i}.  If missing or illegitimate responses exist
    for item \code{i},  the column must also contain an integer
    greater than \code{M_i} that is used to identify such responoses. 
    Alternatively, the column use \code{NA} for this purpose.
    Because missing and illegible responses are normally
    rare, they are given a different and simpler estimation
    procedure for their surprisal values.}
  \item{scoreList}{Either a list of length n, each containing a vector of 
    length \code{M_i} that assigns numeric weights to the options
    for that item.  
    In the special case of multiple choice items where the 
    correct option has weight 1 and all others weight 0, 
    a single integer can identify the correct answer.
    If all the items are of the multiple 
    type, scoreList may be a numeric vector of length n
    containing the right answer indices.  List object
    scoreList is mandatory because these weights define the
    person scores for the surprisal curve estimation process.}
  \item{noption}{A numeric vector of length \code{n}containing the number of 
    choices for each item.  These should not count missing or illegal choices.
    Although this object might seem redundant, it is needed
    for checking the consistencies among other objects and
    as an aid for detecting missing and illegal choices.}
  \item{sumscr_rng}{A numeric vector of length two containing the initial and   
    final values for the interval over which test scores are to be plotted.   
    Default is minimum and maximum sum score.}
  \item{titlestr}{A title string for the data and their analyses.
    Default is NULL.}
  \item{itemlabvec }{A character value containing labels for the items.
    Default is NULL and item position numbers are used.}
  \item{optlabList}{A list vector of length \code{n}, each element \code{i} of 
    which is a character vector of length \code{M_i}.
    Default is NULL, and option numbers are used.}
  \item{nbin}{The number of bins for containing proportions of examinees 
    choosing options.  The default is computed by a function that uses the 
    number of examinees.}
  \item{NumBasis}{The number of spline basis functions used to represent 
    surprisal curves. The default is computed by a function that uses the 
    number of examinees.  
    }
  \item{jitterwrd}{A boolian constant: TRUE implies adding a small random value 
    to each sum score value prior to computing percent rank values.}
  \item{PcntMarkers}{Used in plots of curves to display marker or reference 
    percentage points for abscissa values in plots.}
  \item{verbose}{If TRUE details of calculations are displayed.}
}
\details{
  The score range defined \code{scrrng} should contain all of the sum score 
  values, but can go beyond their boundaries if desired.  For example, 
  it may be that no examinee gets a zero sum score, but for reporting and 
  display purposes using zero as the lower limit seems desirable.
  
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
  \item{chcemat:}{A matrix of response data with N rows and n columns where
    N is number of examinees or respondents and n is number of items.
    Entries in the matrices are the indices of the options chosen.
    Column i of chcemat is expected to contain only the integers 
    \code{1,...,noption}.}  
  \item{optList:}{A list vector containing the numerical score values 
    assigned to the options for this question.}
  \item{key:}{If the data are from a test of the multiple choices type
    where the right answer is scored 1 and the wrong answers 0, this is 
    a numeric vector of length n containing the indices the right answers.  
    Otherwise, it is NULL.}
  \item{Sfd:}{A fd object for the defining the surprisal curves.}
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
  \item{Sdim:}{The total dimension of the surprisal scores.}
  \item{PcntMarkers:}{The marker percentages for plotting: 
    5, 25, 50, 75 and 95.}
  \item{grbg:}{A logical vector of length number of questions. 
    TRUE for an item indicates that a garbage option must be added to the  
    score values, and FALSE indicates that there are no illegal or missing 
    responses and the number of options is equal to number of score values.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{TG_analysis},
        \link{Analyze},
        \link{index_distn},
        \link{index2info},
        \link{index_fun},
        \link{Sbinsmth}}
}
\examples{
  #  Example 1:  Input choice data and key for the short version of the 
  #  SweSAT quantitative multiple choice test with 24 items and 1000 examinees
  #  input the choice data as 1000 strings of length 24
  #  set up index and key data
  chcemat <- Quant_13B_problem_chcemat
  key     <- Quant_13B_problem_key
  # number of examinees and of items
  N <- nrow(chcemat)
  n <- ncol(chcemat)
  # number of options per item and option weights
  noption <- rep(0,n)
  for (i in 1:n) noption[i]  <- 4
  scoreList <- list() # option scores
  for (item in 1:n){
    scorei <- rep(0,noption[item])
    scorei[Quant_13B_problem_key[item]] <- 1
    scoreList[[item]] <- scorei
  }
  # Use the input information to define the 
  # big three list object containing info about the input data
  dataList <- make_dataList(chcemat, scoreList, noption)
  }

