\name{movie}
\alias{movie}
\docType{data}
\title{Movies and Aggressive Affect}
\description{Participants are randomly assigned to one of two groups. The first group watches a violent film, and the other watches a nonviolent film. Afterwards, the aggressive affect is measured, and it is desired to compare three groups, taking gender and degress into account as well.}
\usage{movie}
\format{
  A data frame with 4 variables and 68 observations:
  \describe{
    \item{\code{degree}}{no degree vs. degree}
    \item{\code{gender}}{36 males, 32 females}
    \item{\code{type}}{violend vs. nonviolent}
    \item{\code{aggressive}}{aggressive affect}
  }
}
\details{
Artificial dataset from Wilcox book (p. 316).   
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\examples{
movie
summary(movie)
}
\keyword{datasets}
