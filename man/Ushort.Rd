\name{Ushort}
\alias{Ushort}
\docType{data}
\title{
   Test data for 24 math calculation questions from the SweSAT Quantitative data.
}
\description{
   These data are for a randomly selected subset of 1000 examinees.
}
\usage{
  Ushort
}
\format{
  A .txt file with 1000 rows each containing 24 integers in string format.  The
  integers indicate which answer was chosen for each question by the examinee
  associated with the row.
}
\details{
  The code above inputs the contents of the file using the scan function into
  a single character vector of length 24,000.  The strings are converted
  to integers by the third command, and the integer vector is then
  reformatted into a 1000 by 24 matrix containing choice indices.
  The subsequent analysis of these data is described in detail in the vignette 
  \code{SweSATQuanitativeShort}.
}
\keyword{datasets}
