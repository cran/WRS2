\name{chile}
\alias{chile}
\docType{data}
\title{Chile Heat and Length}
\description{
Originally from pepperjoe.com, this dataset contains the name, length, and heat of chiles. Heat is measured on a scale from 0-11. 
(0-2 ... for sissys, 3-4 ... sort of hot, 5-6 ... fairly hot, 7-8 ... real hot, 9.5-9 ... torrid, 9.5-11 ... nuclear). 
}
\usage{chile}
\format{ 
  A data frame with 3 variables and 85 observations:
  \describe{
    \item{\code{name}}{name of the chile}
    \item{\code{length}}{length in cm}
    \item{\code{heat}}{heat of the chile}
  }
}
\references{
Wright, D. B., & London, K. (2009). Modern Regression Techniques Using R. Sage.
}
\examples{
summary(chile)
}
\keyword{datasets}
