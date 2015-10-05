\name{electric}
\alias{electric}
\docType{data}
\title{The Electric Company}
\description{
These data are based on an educational TV show for children called "The Electric Company". In each of four grades, the classes were 
randomized into treated (TV show) and control groups (no TV show). At the beginning and at the end of the school year, students in all the classes were given a reading test. The average test scores per class were recorded. 
}
\usage{electric}
\format{ 
  A data frame with 5 variables and 192 observations:
  \describe{
    \item{\code{City}}{Fresno and Youngstown}
    \item{\code{Grade}}{grade}
    \item{\code{Pretest}}{reading scores at the beginning of the semester}
    \item{\code{Posttest}}{reading scores at the end of the semester}
    \item{\code{Group}}{treatment vs. control}
  }
}
\references{
Gelman, A., & Hill, J. (2007). Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge University Press: New York, NY. 
}
\examples{
summary(electric)
}
\keyword{datasets}
