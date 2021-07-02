\name{t2way}
\alias{t2way}
\alias{pbad2way}
\alias{mcp2a}
\alias{mcp2atm}
\alias{print.t2way}
\alias{print.mcp}
\alias{model.matrix.mcp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A two-way ANOVA for trimmed means, M-estimators, and medians.
}
\description{
The \code{t2way} function computes a two-way ANOVA for trimmed means with interactions effects. Corresponding post hoc tests are in \code{mcp2atm}. \code{pbad2way} performs a two-way ANOVA using M-estimators for location with \code{mcp2a} for post hoc tests.
}
\usage{
t2way(formula, data, tr = 0.2, ...)
pbad2way(formula, data, est = "mom", nboot = 599, pro.dis = FALSE, ...)
mcp2atm(formula, data, tr = 0.2, ...)
mcp2a(formula, data, est = "mom", nboot = 599, ...)
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
 \item{est}{
Estimate to be used for the group comparisons: either \code{"onestep"} for one-step M-estimator of location using Huber's Psi, \code{"mom"} for the modified one-step (MOM) estimator of location based on Huber's Psi, or \code{"median"}.
}
 \item{nboot}{
number of bootstrap samples.
}
\item{pro.dis}{
If \code{FALSE}, Mahalanobis distances are used; if \code{TRUE} projection distances are computed.
}
\item{...}{
currently ignored.
}
}
\details{
\code{t2way} does not report any degrees of freedom since an adjusted critical value is used.
\code{pbad2way} returns p-values only; if it happens that the variance-covariance matrix in the Mahalanobis distance computation is singular, it is suggested to use the projection distances by setting \code{pro.dis = TRUE}.
}
\value{
The functions \code{t2way} and \code{pbad2way} return an object of class \code{t2way} containing:

\item{Qa}{first main effect}
\item{A.p.value}{p-value first main effect}
\item{Qb}{second main effect}
\item{B.p.value}{p-value second main effect}
\item{Qab}{interaction effect}
\item{AB.p.value}{p-value interaction effect}
\item{call}{function call}
\item{varnames}{variable names}
\item{dim}{design dimensions}

The functions \code{mcp2atm} and \code{mcp2a} return an object of class \code{mcp} containing:

\item{effects}{list with post hoc comparisons for all effects}
\item{contrasts}{design matrix}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{t1way}}, \code{\link{med1way}}, \code{\link{t2way}}
}
\examples{
## 2-way ANOVA on trimmed means
t2way(attractiveness ~ gender*alcohol, data = goggles)

## post hoc tests
mcp2atm(attractiveness ~ gender*alcohol, data = goggles)

## 2-way ANOVA on MOM estimator
pbad2way(attractiveness ~ gender*alcohol, data = goggles)

## post hoc tests
mcp2a(attractiveness ~ gender*alcohol, data = goggles)

## 2-way ANOVA on medians
pbad2way(attractiveness ~ gender*alcohol, data = goggles, est = "median")

## post hoc tests
mcp2a(attractiveness ~ gender*alcohol, data = goggles, est = "median")

## extract design matrix
model.matrix(mcp2a(attractiveness ~ gender*alcohol, data = goggles, est = "median"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

