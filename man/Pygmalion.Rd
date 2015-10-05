\name{Pygmalion}
\alias{Pygmalion}
\docType{data}
\title{Pygmalion Data}
\description{
The Pygmalion effect is the phenomenon where higher expectations lead to an increase in performance. For instance, when teachers expect students to do well and show intellectual growth, they do; when teachers do not have such expectations, performance and growth are not so encouraged and may in fact be discouraged in a variety of ways. This dataset contains reasoning IQ scores of children. For the experimental group, positive expectancies 
had been suggested to teachers after the pretest. For the experimental group, no expectancies had been suggested after the pretest. 
For both groups we have reasoning IQ posttest scores. The dataset is taken from Elashoff and Snow (1970). 
}
\usage{Pygmalion}
\format{
  A data frame with 3 variables and 114 observations:
  \describe{
    \item{\code{Pretest}}{pretest score}
    \item{\code{Posttest}}{posttest score}
    \item{\code{Group}}{treatment vs. control}
  }
}
\references{
Elashoff, J. D., & Snow, R. E. (1970). A case study in statistical inference: Reconsideration of the Rosenthal-Jacobson data on teacher 
expectancy. Technical Report No. 15, School of Education, Stanford University.

Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\examples{
summary(Pygmalion)
}
\keyword{datasets}
