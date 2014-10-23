\name{picture}
\alias{picture}
\docType{data}
\title{Profile Pictures}
\description{This dataset examines how the profile pictures on social network platforms affect the number of friend requests when females are in a relationship. The relationship status is a between-subject variable (part of the participants did set their status to relationship). For the first 3 weeks the subjects had a picture of their own in their profiles. For the following 3 weeks they posted a picture with a man.}
\usage{picture}
\format{
  A data frame with 4 variables and 40 observations:
  \describe{
    \item{\code{case}}{subject id}
    \item{\code{relationship_status}}{Relationship status on social network platform}
    \item{\code{couple}}{amount of friend requests when profile picture as couple}
    \item{\code{alone}}{amount of friend requests when profile picture as single}
  }
}
\details{
Dataset from Field et al. book (p. 644).   
}
\references{
Field, A., Miles, J, \& Field, Z. (2012). Discovering Statistics Using R. Sage.
}
\examples{
picture
summary(picture)
}
\keyword{datasets}
