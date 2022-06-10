\name{invisibility}
\alias{invisibility}
\docType{data}
\title{Cloaks of Invisibility}
\description{We are interested in the effect that wearing a cloak of invisibility has on people's tendency to mischief. 80 participants were placed in an enclosed community. Hidden cameras recorded mischievous acts. It was recorded how many mischievous acts were conducted in the first 3 weeks (mischief1). After 3 weeks 34 participants were told that the cameras were switched off so that no one would be able to see what they're getting up to. The remaining 46 subjects were given a cloak of invisibility. These people were told not to tell anyone else about their cloak and they could wear it whenever they liked. The number of mischievous acts were recorded over the next 3 weeks (mischief2).}
\usage{invisibility}
\format{
  A data frame with 3 variables and 80 observations:
  \describe{
    \item{\code{cloak}}{factor with 34 subjects in the no cloak condition, 46 in the cloak condition}
    \item{\code{mischief1}}{number of mischievous acts during the first 3 weeks}
    \item{\code{mischief2}}{number of mischievous acts during the second 3 weeks}
  }
}
\details{
Fictional dataset from Field et al. book (p. 485).   
}
\references{
Field, A., Miles, J., & Field, Z. (2012). Discovering Statistics Using R. Sage.
}
\examples{
invisibility
summary(invisibility)
}
\keyword{datasets}
