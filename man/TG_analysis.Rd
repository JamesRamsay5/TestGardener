\name{TG_analysis}
\alias{TG_analysis}
\title{Statistics for Multiple choice Tests, Rating Scales and Other Choice Data)
}
\description{
  Given an choice ingter-valued index matrix and a vector of numbers 
  of item options, the function cycles through a set of iterations 
  involving surprisal curve estimation followed by test taker index 
  values.
}
\usage{
  TG_analysis(chcemat, scoreList, noption, NumBasis=7, ncycle=10,  
              titlestr=NULL, itemlabvec=NULL, optlabList=NULL,
              nbin=nbinDefault(N), sumscr_rng=NULL,
              jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
              cyclechoice=ncycle, itdisp=FALSE, verbose=FALSE)
}
\arguments{
  \item{chcemat}{An N by n matrix.  Column i must contain the integers 
  from 1 to M_i, where M_i is the number of options
  for item i.  If missing or illegitimate responses exist
  for item i,  the column must also contain an integer
  greater than M_i that is used to identify such responoses. 
  Alternatively, the column use NA for this purpose.
  Because missing and illegible responses are normally
  rare, they are given a different and simpler estimation
  procedure for their surprisal values.
  U is mandatory.}
  \item{scoreList}{Either a list of length n, each containing a vector of 
  length M_i that assigns numeric weights to the options
  for that item.  
  In the special case of multiple choice items where the 
  correct option has weight 1 and all others weight 0, 
  a single integer can identify the correct answer.
  If all the items are of the multiple 
  type, scoreList may be a numeric vector of length n
  containing the right answer indices.  List object
  scoreList is mandatory because these weights define the
  person scores for the surprisal curve estimation process.}
  \item{noption}{A numeric vector of length \code{n} containing the 
  number of options for each item.}
  \item{NumBasis}{The number of spline basis functions to use for 
  surprisal values.  Defaults to 7.}
  \item{ncycle}{The number of cycles in the analysis.  Defaults to 10.}
  \item{titlestr}{A title string for the data and their analyses.
  Default is NULL.}
  \item{itemlabvec}{A character value containing labels for the items.
  Default is NULL and item position numbers are used.}
  \item{optlabList}{A list vector of length n, each element i of which is a
  character vector of length M_i.
  Default is NULL, and option numbers are used.}
  \item{nbin}{The number of bins containing proportions of choices.}
  \item{sumscr_rng}{A vector of length 2 indicating the initial and final
  sum score values.  Default is NULL the whole sum score
  is used.}
  \item{jitterwrd}{A logical object indicating whether a small jittering
  perturbation should be used to break up ties.  
  Defaults to TRUE.}
  \item{PcntMarkers}{A vector of percentages inside of [0,100] that appear
  in plots.  Defaults to c(5, 25, 50, 75, 95).
  Extra displays are provided.  Defaults to FALSE.}
  \item{itdisp}{Display results for function \code{theta_fun}.}
  \item{verbose}{Extra displays are provided.  Defaults to FALSE.}
  \item{cyclechoice}{A number within 1 to 10 indicating which cycle will be 
  used to represent the TestGardener results.  
  Defaults to ncycle.}
}
\details{
  This function in package TestGardener processes at a mininum two objects:  
  (1) A matrix \code{chcemat} that contains indices of choices made in a 
  sequence of choice situations (its number columns \code{n}) by a set of 
  persons making the choices (its number of rows \code{N}); and
  (2) A list vector \code{scoreList} of length \code{n} containing numerical 
  weights or scores for each choice available with in each of \code{n} choice
  situations (referred to as \code{items}).
  
  The function returns three large lists containing objects that can be used
  to assess:
  (1) the probability that a choice will be made, and
  (2) the quantity of information, called \code{surprisal}, that the choice 
  made reveals about the performance or experience of the person making the 
  choice.
}
\value{
Four list objects, each containing objects that are required for various
displays, tables and other results:
  \item{parmList}{A list object containing objects useful for displaying
    results that involve the score index cotninuum:}
    \itemize{
      \item{SfdList: A list object of length n, each containing  
      objects for an item for displaying that item's surprisal curves
      as defined by the score index values after the analysis.
      See the help page for function \code{Analyze} for a
      description of these objects.}
      \item{Qvec: A vector containing the positions on the score index
      continuum of the marker percentages defined in the arguments of
      function \code{make_dataList()}.}
      \item{binctr: A vector of length \code{nbin} containing the positions 
      on the score index continuum of the bin centres.}
      \item{indexScore: A vector of length \code{N} containing 
      the positions on the score index contiuum of each person.}
      \item{infoSurp: The length of the test or scale information 
      continuum in M-bits.}
    }
  \item{infoList}{A list object containing objects useful for displaying
    results that involve the scale information cotninuum:}
    \itemize{
      \item{infofine: A fine mesh of 101 values that is used to plot
      the scale information comntinuum.}
      \item{scopevec: A vector of length \code{N} containing 
      the positions on the scale information contiuum of each person.}
      \item{Qinfovec: A vector containing the positions on the scale
      information continuum of the marker percentages defined in the 
      arguments of function \code{make_dataList()}.}
      \item{infobinctr: A vector of length \code{nbin} containing the positions 
      on the scale information continuum of the bin centres.}
    }
  \item{HALsave}{A \code{ncycle} by 2 matrix containing in the first column
  the data fit values and in the second column the infoSurp or arclength
  values.}
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
        \link{Analyze},
        \link{index_distn},
        \link{index2info},
        \link{index_fun},
        \link{Sbinsmth}}
}
\examples{
\dontrun{
  #  Example 1:  Input choice data and key for the short version of the 
  #  SweSAT quantitative multiple choice test with 24 items and 1000 examinees
  #  input the choice indices in the 1000 by 24 choice index matrix
  chcemat <- Quant_13B_problem_chcemat
  #  set up the key data
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
    scorei[key[item]] <- 1
    scoreList[[item]] <- scorei
  }
  #  Analyze the data and return the big three list objects:
  #  dataList (input info), parmList (score index info) and
  #  infoList (information or scope info)
  TGresult <- TG_analysis(chcemat, scoreList, noption, 
                          NumBasis=4, ncycle=10, verbose=TRUE)}
}

