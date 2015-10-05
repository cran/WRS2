\name{yuen}
\alias{yuen}
\alias{yuenbt}
\alias{yuenbt}
\alias{pb2gen}
\alias{print.yuen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Independent samples t-tests on robust location measures.
}
\description{
The function \code{yuen} performs Yuen's test for trimmed means, \code{yuenbt} is a bootstrap version of it. The \code{pb2gen} function performs a t-test based on various robust estimators. 
}
\usage{
yuen(formula, data, tr = 0.2)
yuenbt(formula, data, tr = 0.2, nboot = 599)
pb2gen(formula, data, est = "mom", nboot = 599)
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
\item{est}{
Estimate to be used for the group comparisons: either \code{"onestep"} for one-step M-estimator of location using Huber's Psi, \code{"mom"} for the modified one-step (MOM) estimator of location based on Huber's Psi, or \code{"median"}, \code{"mean"}. 
}
}


\value{
Returns objects of classes \code{"yuen"} or \code{"pb2"} containing:

\item{test}{value of the test statistic (t-statistic)}
\item{p.value}{p-value}
\item{conf.int}{confidence interval}
\item{df}{degress of freedom}
\item{diff}{trimmed mean difference}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.

Yuen, K. K. (1974). The two sample trimmed t for unequal population variances. Biometrika, 61, 165-170.
}
\seealso{
\code{\link{t1way}},\code{\link{t1waybt}}
}
\examples{
## Yuen's test
yuen(Anxiety ~ Group, data = spider)

## Bootstrap version of Yuen's test (symmetric CIs)
yuenbt(Anxiety ~ Group, data = spider)

## Using an M-estimator
pb2gen(Anxiety ~ Group, data = spider, est = "mom")
pb2gen(Anxiety ~ Group, data = spider, est = "mean")
pb2gen(Anxiety ~ Group, data = spider, est = "median")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

