\name{diet}
\alias{diet}
\docType{data}
\title{Diet and Weight Loss}
\description{
Weight loss is studied for three different types of diets.
}
\usage{bush}
\format{
  A data frame with 7 variables and 78 observations:
  \describe{
    \item{\code{gender}}{gender}
    \item{\code{age}}{age}
    \item{\code{height}}{body height}
    \item{\code{diet.type}}{three types of diet}
    \item{\code{initial.weight}}{initial weight before diet}
    \item{\code{final.weight}}{final weight after diet}
    \item{\code{weight.loss}}{weight loss}
  }
}
\details{
Original data taken from \url{https://www.sheffield.ac.uk/polopoly_fs/1.570199!/file/stcp-Rdataset-Diet.csv}, slightly modified.
}
\references{
Couturier, D. L., Nicholls, R., and Fernandes, M. (2018). ANOVA with R: analysis of the diet dataset. Retrieved from 
\url{https://bioinformatics-core-shared-training.github.io/linear-models-r/ANOVA.html#section_1:_importation_and_descriptive_analysis}
}

\examples{
str(diet)
}
\keyword{datasets}
