% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyzer_alpha_0.06.R
\name{sign_detect}
\alias{sign_detect}
\title{Sign detection}
\usage{
sign_detect(dat)
}
\arguments{
\item{dat}{the vector to be used}
}
\value{
vector with the signs of each element of the original matrix
}
\description{
This function takes as input a vector and returns the sign of each element
}
\examples{
sign_detect(c(1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1))

}
