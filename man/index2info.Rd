\name{index2info}
\alias{index2info}
\title{
  Compute results using arc length or information as the abscissa.
}
\description{The one-dimensional psychometric model defines a space curve within 
  the vector space defined by the total collection of option surprisal curves.  
  This curve is a valuable resource since positions along the curve are defined 
  in bits and positions on the curve are subject to the same strict properties 
  that apply to physical measurements.

  Function \code{index2info} is required to convert
  objects defined over the score index continuum \code{c(0,100)} to
  the same objects over the arc length continuum \code{c(0,infoSurp)}, and also
  vice versa.  Since the arc length or information continuum is along
  a space curve that is invariant under strictly monotone transformations
  of the score index \code{index}, and is also a metric, it is an ideal 
  choice for the abscissa in all plots.}
\usage{
   index2info(index, Qvec, SfdList, binctr, itemindex=1:n, plotrng=c(0,100), 
                shortwrd)
}
\arguments{
  \item{index}{A vector of score index, test score, or arc length values, one 
    for each examinee or respondent.}
  \item{Qvec}{A vector of locations of the five marker percentages.}
  \item{SfdList}{A numbered list object produced by a TestGardener analysis of 
  a test.  Its length is equal to the number of items in the test or questions
  in the scale.  Each member of \code{SfdList} is a named list containing 
  information computed during the analysis.}
  \item{binctr}{A vector of locations of the bin centers.}
  \item{itemindex}{A vector containing the indices of the items to be used.}
  \item{plotrng}{A vector of length 2 containing the starting score index
    and end score index values of the range to be plotted.}
  \item{shortwrd}{If TRUE only vectors infoSurp and infoSurpvec are returned in 
    order to speed up the computation within cycles in function \code{Analyze()} 
    where only these objects are required.  The default is FALSE.}
}
\value{A named list object containing these results of the analysis:
  \item{infoSurp}{The length of the test information or scale curve.}
  \item{infoSurpvec}{Positions on the test information or scale curve 
    corresponding to a fine mesh of score index values (typically 101 values 
    between 0 and 100).} 
  \item{infoSurpfd}{Functional data object representing the relation between
     the score index abscissa and the infoSurp or information ordinate.}
  \item{scopevec}{A vector of positions on the test information or scale curve
    corresponding to the input score index values in argument \code{index}.} 
  \item{Qvec_al}{Values in arc length of the five marker percentages.}
  \item{binctr_al}{Values in arc length of the bin centers.}
  \item{Sfd.info}{A functional data object representing the relation between
    the infoSurp or information abscissa and the score index ordinate.}
  \item{Sdim.index}{The dimension of the overspace, which equal to sum of the 
    number of options in the items specified in \code{itemindex}.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Analyze}}
}
\examples{
  #  Example 1.  Display the scope or information curve for the 
  #  short SweSAT multiple choice test with 24 items and 1000 examinees.
  #  The scope curve is constructed using the complete analysis cycles.
  #  Set up the required arguments using the converged parmList object.
  indfine     <- seq(0,100,len=101)
  index       <- Quant_13B_problem_parmList$index
  Qvec        <- Quant_13B_problem_parmList$Qvec
  SfdList     <- Quant_13B_problem_parmList$SfdList
  binctr      <- Quant_13B_problem_parmList$binctr
  #  Carry out the construction of the information results.
  infoList    <- index2info(index, Qvec, SfdList, binctr)
  # Plot the shape of the information curve
  oldpar <- par(no.readonly=TRUE)
  Scope_plot(infoList$infoSurp, infoList$infoSurpvec)
  par(oldpar)
}
