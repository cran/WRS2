\name{hangover}
\alias{hangover}
\docType{data}
\title{Hangover Symptoms}
\description{
In a study on the effect of consuming alcohol, hangover symptoms were measured for two independent groups, with each subject consuming alcohol and being measured on three different occasions. One group consisted of sons of alcoholics and the other was a control group.
}
\usage{hangover}
\format{ 
  A data frame with 4 variables and 120 observations:
  \describe{
    \item{\code{symptoms}}{number of hangover symptoms}
    \item{\code{group}}{son of alcoholic vs. control}
    \item{\code{time}}{measurement occasion}
    \item{\code{id}}{subject ID}
  }
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\examples{
summary(hangover)
}
\keyword{datasets}
