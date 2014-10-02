\name{med1way}
\alias{med1way}
\alias{print.med1way}
\title{
Heteroscedastic one-way ANOVA for medians.
}
\description{
This function computes a one-way ANOVA for the medians. Homoscedasticity assumption not required. 
}
\usage{
med1way(formula, data, alpha = 0.05, crit = NA, iter = 1000, SEED = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
  \item{alpha}{
significance level. 
}
  \item{crit}{
if \code{NA} an appropriate critical value is determined for the alpha value specified and a p-value is returned. If a value for \code{crit} is specified, it is used as the critical value, but no p-value is reported. Specifying a value for \code{crit} reduces execution time.
}
  \item{iter}{
number of iterations.
}
  \item{SEED}{
whether a seed for the random numbers should be set.
}
}
\value{
Returns an object of class \code{med1way} containing:

\item{test}{value of the test statistic}
\item{crit.val}{critical value}
\item{p.value}{p-value}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{t1way}},\code{\link{t1waybt}}
}
\examples{
med1way(libido ~ dose, data = viagra)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
