\name{WineTasting}
\alias{WineTasting}
\docType{data}
\title{Wine Tasting}
\description{In this hypothetical dataset we have three types of wine (A, B and C). We asked 22 friends to taste each of the three wines (in a blind fold fashion), and then to give a grade of 1 to 7. We asked them to rate the wines 5 times each, and then averaged their results to give a number for a persons preference for each wine. 
}
\usage{WineTasting}
\format{
  A data frame with 3 variables and 66 observations:
  \describe{
    \item{\code{Taste}}{Taste Rating}
    \item{\code{Wine}}{Wine (A, B, C)}
    \item{\code{Taster}}{Taster (index)}
  }
}
\examples{
WineTasting
summary(WineTasting)
}
\keyword{datasets}
