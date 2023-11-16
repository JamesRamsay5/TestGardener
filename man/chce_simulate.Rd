\name{chce_simulate}
\alias{chce_simulate}
\alias{chce_simulate}
\title{
  Simulate a test or scale data matrix.
}
\description{
  Used in \code{dataSimulation}, this function sets up an N by n matrix of 
  index values that specify the index of the option chosen by an examinee or
  respondent for a specific question.
}
\usage{
  chce_simulate(index.pop, SfdList)
}
\arguments{
  \item{index.pop}{A vector containing population score index values at which 
    data are to be simulated.}
  \item{SfdList}{A numbered list object produced by a TestGardener analysis of 
  a test.  Its length is equal to the number of items in the test or questions 
  in the scale.  Each member of \code{SfdList} is a named list containing 
  information computed during the analysis.}
}
\details{For each question and each examinee a vector of random multinomial 
integer values is generated using the probability transforms of the surprisal 
curves and the examinee's score index value.}
\value{An N by n matrix of integer index values.}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
