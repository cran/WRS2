\name{WineTasting}
\alias{WineTasting}
\docType{data}
\title{Wine Tasting}
\description{In this hypothetical dataset we have three types of wine (A, B and C), and we would like to know which one is the best one (in a scale of 1 to 7). We asked 22 friends to taste each of the three wines (in a blind fold fashion), and then to give a grade of 1 till 7 (we asked them to rate the wines 5 times each, and then averaged their results to give a number for a persons preference for each wine. This number which is now an average of several numbers.
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
