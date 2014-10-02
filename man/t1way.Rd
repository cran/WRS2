\name{t1way}
\alias{t1way}
\alias{print.t1way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A heteroscedastic one-way ANOVA for trimmed means.
}
\description{
This function computes a one-way ANOVA for the medians. Homoscedasticity assumption not required. It uses a generalization of Welch's method. 
}
\usage{
t1way(formula, data, tr = 0.2)
}
%- maybe also 'usage' for other objects documented here.
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
\value{
Returns an object of class \code{t1way} containing:

\item{test}{value of the test statistic (F-statistic)}
\item{df1}{degrees of freedom}
\item{df2}{degrees of freedom}
\item{p.value}{p-value}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{med1way}},\code{\link{t1waybt}}
}
\examples{
t1way(libido ~ dose, data = viagra)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

