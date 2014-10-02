\name{goggles}
\alias{goggles}
\docType{data}
\title{Beer Goggles Effect}
\description{This dataset is about the effects of alcohol on mate selection in night-clubs. The hypothesis is that after alcohol had been consumed, subjective perceptions of physical attractiveness would become more inaccurate (beer-goggles effect).}
\usage{goggles}
\format{
  A data frame with 3 variables and 48 observations:
  \describe{
    \item{\code{gender}}{24 male, 24 female students}
    \item{\code{alcohol}}{amount of alcohol consumed}
    \item{\code{attractiveness}}{attractiveness rating (0-100)}
  }
}
\details{
Dataset from Field et al. book (p. 501).   
}
\references{
Field, A., Miles, J, \& Field, Z. (2012). Discovering Statistics Using R. Sage.
}
\examples{
goggles
summary(goggles)
}
\keyword{datasets}
