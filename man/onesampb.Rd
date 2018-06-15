\name{onesampb}
\alias{onesampb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{One-sample percentile bootstrap}
\description{
Tests whether a robust location measure (median, Huber Psi) differs from a null value and reports a 95\% confidence interval based on percentile bootstrap. 
}
\usage{
onesampb(x, est = "onestep", nboot = 2000, nv = 0)
}
\arguments{
  \item{x}{
a numeric vector.
}
 \item{est}{
robust estimator to be used (\code{"onestep"}, \code{"mom"}, or \code{"median"}).
}
  \item{nboot}{
number of bootstrap samples.
}
 \item{nv}{
value under H0. 
 }
}
\value{
Returns an object of class \code{"onesampb"} containing:

\item{ci}{95\% confidence interval}
\item{estimate}{robust location sample estimate}
\item{p.value}{p-value}
\item{n}{number of effective observations}
\item{call}{function call}
}
\references{
Wilcox, R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.). Elsevier.
}
\seealso{
\code{\link{t1way}}
}
\examples{
set.seed(123)
x <- rnorm(30)
onesampb(x, nboot = 100)    ## H0: Psi = 0

set.seed(123)
x <- rlnorm(30)
onesampb(x, est = "median", nv = 1)   ## H0: median = 1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

