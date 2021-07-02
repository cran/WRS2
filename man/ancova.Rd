\name{ancova}
\alias{ancova}
\alias{ancboot}
\alias{print.ancova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust ANCOVA
}
\description{
This function computes robust ANCOVA for 2 independent groups and one covariate. It compares trimmed means. No parametric assumption (e.g. homogeneity) is made about the form of the regression lines. A running interval smoother is used. A bootstrap version which computes
confidence intervals using a percentile t-bootstrap is provided as well.
}
\usage{
ancova(formula, data, tr = 0.2, fr1 = 1, fr2 = 1, pts = NA, ...)

ancboot(formula, data, tr = 0.2, nboot = 599, fr1 = 1, fr2 = 1, pts = NA, ...)
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
 \item{fr1}{
values of the span for the first group (1 means unspecified)
}
 \item{fr2}{
values of the span for the second group (1 means unspecified)
}
 \item{pts}{
can be used to specify the design points where the regression lines are to be compared; if \code{NA} design points are chosen.
}
\item{nboot}{
number of bootstrap samples
}
\item{...}{
currently ignored.
}
}
\value{
Returns an object of class \code{ancova} containing:

\item{evalpts}{covariate values (including points close to these values) where the test statistic is evaluated}
\item{n1}{number of subjects at evaluation point (first group)}
\item{n2}{number of subjects at evaluation point (first group)}
\item{trDiff}{trimmed mean differences}
\item{se}{standard errors for trimmed mean differences}
\item{ci.low}{lower confidence limit for trimmed mean differences}
\item{ci.hi}{upper confidence limit for trimmed mean differences}
\item{test}{values of the test statistic}
\item{crit.vals}{critical values}
\item{p.vals}{p-values}
\item{fitted.values}{fitted values from interval smoothing}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{t2way}}
}
\examples{
head(invisibility)
ancova(mischief2 ~ cloak + mischief1, data = invisibility)

## specifying covariate evaluation points
ancova(mischief2 ~ cloak + mischief1, data = invisibility, pts = c(3, 4, 8, 1))

## bootstrap version
ancboot(mischief2 ~ cloak + mischief1, data = invisibility)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

