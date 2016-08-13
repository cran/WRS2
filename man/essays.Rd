\name{essays}
\alias{essays}
\docType{data}
\title{Academic Writing Data}
\description{This study looked at the effects of two forms of written corrective feedback on lexico-grammatical accuracy in the academic writing of English as a foreign language university students. It had a 3x4 within-by-between design with three groups (two treatment and one control) measured over four occasions (pretest/treatment, treatment, posttest, delayed posttest).
}
\usage{bush}
\format{
  A data frame with 4 variables and 120 observations:
  \describe{
    \item{\code{id}}{participant ID}
    \item{\code{group}}{control, direct, indirect}
    \item{\code{essay}}{four measurement occasions}
    \item{\code{errorRatio}}{error ratio}
  }
}

\references{
McGrath, D. (2016). The Effects of Comprehensive Direct and Indirect Written Corrective Feedback on Accuracy in English as a Foreign Language Students' Writing (Unpublished master's thesis). Macquarie University, Sydney, Australia.
}
\examples{
head(essays)
summary(essays)
}
\keyword{datasets}
