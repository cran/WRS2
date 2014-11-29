\name{t1way}
\alias{t1way}
\alias{lincon}
\alias{print.t1way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A heteroscedastic one-way ANOVA for trimmed means.
}
\description{
The \code{t1way} function computes a one-way ANOVA for the medians. Homoscedasticity assumption not required. It uses a generalization of Welch's method. Corresponding post hoc tests can be performed using \code{lincon}.
}
\usage{
t1way(formula, data, tr = 0.2)
lincon(formula, data, tr = 0.2)
}
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
  \item{tr}{
trim level for the mean.
}
}
\details{
In the post hoc computations, confidence intervals are adjusted to control FWE, but p-values are not adjusted to control FWE. 
}
\value{
\code{t1way} returns an object of class \code{"t1way"} containing:

\item{test}{value of the test statistic (F-statistic)}
\item{df1}{degrees of freedom}
\item{df2}{degrees of freedom}
\item{p.value}{p-value}
\item{call}{function call}

\code{lincon} returns an object of class \code{"mcp1"} containing:

\item{comp}{inference for all pairwise comparisons}
\item{fnames}{names of the factor levels}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{med1way}},\code{\link{t1waybt}}
}
\examples{
t1way(libido ~ dose, data = viagra)

## post hoc tests
lincon(libido ~ dose, data = viagra)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

