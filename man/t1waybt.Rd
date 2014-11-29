\name{t1waybt}
\alias{t1waybt}
\alias{mcppb20}
\alias{print.t1waybt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap version of the heteroscedastic one-way ANOVA for trimmed means. 
}
\description{
Test the hypothesis of equal trimmed means using a percentile t bootstrap method. Corresponding post hoc tests are provided in \code{mcppb20}.
}
\usage{
t1waybt(formula, data, tr = 0.2, nboot = 599)
mcppb20(formula, data, tr = 0.2, nboot = 599, crit = NA)
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
  \item{nboot}{
number of bootstrap samples. 
}
  \item{crit}{
critical significance value.
}
}

\value{
Returns an object of class \code{t1waybt} containing:

\item{test}{value of the test statistic}
\item{p.value}{p-value}
\item{Var.Explained}{explained amount of variance}
\item{Effect.Size}{effect size}
\item{nboot.eff}{effective number of bootstrap samples}
\item{call}{function call}

\code{mcppb20} returns an object of class \code{"mcp1"} containing:

\item{comp}{inference for all pairwise comparisons}
\item{fnames}{names of the factor levels}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{t1way}},\code{\link{med1way}}
}
\examples{
t1waybt(libido ~ dose, data = viagra)

## post hoc
mcppb20(libido ~ dose, data = viagra)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% __ONLY ONE__ keyword per line
