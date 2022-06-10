\name{bush}
\alias{bush}
\docType{data}
\title{Bushtucker Foods}
\description{In the TV show "I'm a celebrity, get me out of here" the celebrities had to eat things like stick insects, fish eyes, etc.
This dataset records the time taken to retch when eating these things.
}
\usage{bush}
\format{
  A data frame with 5 variables and 8 observations:
  \describe{
    \item{\code{participant}}{participant ID}
    \item{\code{stick_insect}}{time taken to retch when eating a stick insect}
    \item{\code{kangaroo_testicle}}{time taken to retch when eating a kangaroo testicle}
    \item{\code{fish_eye}}{time taken to retch when eating a fish eye}
    \item{\code{witchetty_grub}}{time taken to retch when eating a witchetty grub}
  }
}
\details{
Dataset from Field et al. book (p. 557).   
}
\references{
Field, A., Miles, J., & Field, Z. (2012). Discovering Statistics Using R. Sage.
}
\examples{
bush
summary(bush)
}
\keyword{datasets}
