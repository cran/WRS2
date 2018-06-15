\name{med1way}
\alias{med1way}
\alias{print.med1way}
\title{
Heteroscedastic one-way ANOVA for medians.
}
\description{
This function computes a one-way ANOVA for the medians. Homoscedasticity assumption not required. There shouldn't be too many ties.
}
\usage{
med1way(formula, data, iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
\item{iter}{
number of iterations to determine critical value.
}
}
\details{
Evaluating the test statistic using the df proposed in the literature can result in the actual level being less than the nominal level, (i.e., around 0.02-0.025 when testing at the 0.05 level and the sample size is small). A better strategy is to simulate the critical 
value and computing the p-value accordingly, as implemented in this function.  
}
\value{
Returns an object of class \code{med1way} containing:
\item{test}{F-value of the test statistic}
\item{crit.val}{critical value}
\item{p.value}{p-value}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{t1way}},\code{\link{t1waybt}}, \code{\link{Qanova}}
}
\examples{
set.seed(123)
med1way(libido ~ dose, data = viagra, iter = 3000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
