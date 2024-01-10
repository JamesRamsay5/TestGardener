\name{ICC}
\alias{ICC}
\title{Plotting probability and surprisal curves for an item}
\description{This is an S3 object that contains information essential plotting 
probability and surprisal curves for a single multiple choice or rating
question.  Bin probabilities and surprisal values can also be plotted.}
\usage{ICC(x, M, Sfd, Zmat, Pbin, Sbin, Pmatfine, Smatfine, DSmatfine, D2Smatfine,
                PStdErr, SStdErr, ItemArcLen, itemStr=NULL, optStr=NULL)}
\arguments{
  \item{x}{An item number.}
  \item{M}{The number of options for this item, including an option for missing or 
  illegal values if required.}
  \item{Sfd}{A functional surprisal curve object defined by \code{K} B-spline 
  basis functions and a \code{K} by \code{M-1} matrix of coefficients.}
  \item{Zmat}{An \code{M} by \code{M-1} matrix satisfying the conditions
  \code{t(Zmat) Zmat} = \code{I} and columns sum to zero.}
  \item{Pbin}{A \code{nbin} by \code{M} matrix of probabilities that a given
  bin is chosen by a test taker.}
  \item{Sbin}{A \code{nbin} by \code{M} matrix of surprisal values for the 
  probabilities in \code{Pbin.}}
  \item{Pmatfine}{A 101 by \code{M} matrix of probability curve values over
  equally-spaced score index values spanning the interval [0,100].}
  \item{Smatfine}{A 101 by \code{M} matrix of surprisal curve values 
  corresponding to the probability values in \code{Pmatfine}.}
  \item{DSmatfine}{A 101 by \code{M} matrix of first derivative values
  with respect to score index values for the surprisal values.}
  \item{D2Smatfine}{A 101 by \code{M} matrix of second derivative values.}
  \item{PStdErr}{A 101 by \code{M} matrix of standard error estimates for the 
  probability curve values.}
  \item{SStdErr}{A 101 by \code{M} matrix of standard error estimates for the 
  surprisal curve values.}
  \item{ItemArcLen}{The scope or arc length of the item curve.}
  \item{itemStr}{A string that is the name of the item.}
  \item{optStr}{A character vector containing labels for the item options.}
}

\details{The name ICC for this object is an acronym for the term
"item characteristic curve" widely used in the psychometric commuunity.

Function ICC is set up after the initialization process in function 
\code{make_dataList()} has created the members of \code{dataList}.  
Within this list is object \code{SfdList}, which cintains a functional data
object \code{Sfd} for each item.  Both the intial coefficient matrices and
the subsequent estimates of them are available from \code{Sfd$coefs}, and
therefore are available in the ICC object.  These coefficient matrices are
\code{K} by \code{M-1} where \code{K} is the number of basis functions and 
\code{M} is the number of options for asn item.}

\value{The values returned are simply those in the argument list. The S3 
ICC object checks each of these and makes available the S3 commands or methods
\code{str}, \code{print} and \code{plot} that apply the corresponding 
\code{ICC} versions of these opterations.}

\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}

\author{Juan Li and James Ramsay}
