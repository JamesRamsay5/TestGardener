\name{keyshsort}
\alias{keyshort}
\docType{data}
\title{
  Option information for the short form of the SweSAT Quantitative test.
}
\description{
  The file has two rows, each with 24 integers in string format.
  The first row contains the indices of the right answers among the options
  for the 24 questions.  The second row contains the number of options for
  each question.
}
\usage{
  keyshort
}
\details{
  The code above inputs the contents of the file using the scan function into
  a single character vector of length 48.  The strings are converted
  to integers by the second command, and the integer vector is then
  reformatted into a 2 by 24 matrix.  
  The subsequent analysis of these data is described in detail in the vignette 
  \code{SweSATQuanitativeShort}.
}
\keyword{datasets}
