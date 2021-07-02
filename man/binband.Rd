\name{binband}
\alias{binband}
\alias{discANOVA}
\alias{discmcp}
\alias{discstep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparison of discrete distributions
}
\description{
\code{binband} compares two independent variables in terms of their probability function. \code{discANOVA} Tests the global hypothesis that for two or more independent groups, the corresponding discrete distributions are identical. That is, test the hypothesis that independent groups have identical multinomial distributions. \code{discmcp} provides multiple comparisons for J independent groups having discrete distributions. \code{discstep} implements the step-down multiple comparison procedure for comparing J independent discrete random variables.
}
\usage{
binband(x, y, KMS = FALSE, alpha = 0.05, ADJ.P = FALSE, ...)
discANOVA(formula, data, nboot = 500, ...)
discmcp(formula, data, alpha = 0.05, nboot = 500, ...)
discstep(formula, data, nboot = 500, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
an numeric vector of data values for group 1.
}
 \item{y}{
an numeric vector of data values for group 2.
}
 \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
\item{nboot}{
number of bootstrap samples.
}
\item{alpha}{
alpha level.
}
\item{KMS}{
whether the Kulinskaya-Morgenthaler-Staudte method for comparing binomials should be used.
}
\item{ADJ.P}{
whether the critical p-value should be adjusted to control FWE when the sample size is small (<50)
}
\item{...}{
currently ignored.
}
}
\value{
\code{discANOVA} returns an object of class \code{"med1way"} containing:

\item{test}{value of the test statistic}
\item{crit.val}{critical value}
\item{p.value}{p-value}
\item{call}{function call}

The remaining functions return an object of class \code{"robtab"} containing:

\item{partable}{parameter table}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.

Kulinskaya, E., Morgenthaler, S. and Staudte, R. (2010). Variance stabilizing the difference of two binomial proportions. American Statistician, 64, p. 350-356.
}
\seealso{
\code{\link{t1way}}, \code{\link{Qanova}}
}
\examples{
## Consider a study aimed at comparing two methods for reducing shoulder pain after surgery.
## For the first method, the shoulder pain measures are:
x1 <- c(2, 4, 4, 2, 2, 2, 4, 3, 2, 4, 2, 3, 2, 4, 3, 2, 2, 3, 5, 5, 2, 2)
## and for the second method they are:
x2 <- c(5, 1, 4, 4, 2, 3, 3, 1, 1, 1, 1, 2, 2, 1, 1, 5, 3, 5)

fit1 <- binband(x1, x2)
fit1

fit2 <- binband(x1, x2, KMS = TRUE, alpha = 0.01)
fit2

## More than two groups:
discANOVA(libido ~ dose, viagra, nboot = 200)

## Multiple comparisons:
discmcp(libido ~ dose, viagra)

discstep(libido ~ dose, viagra)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

