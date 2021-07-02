\name{twocor}
\alias{twocor}
\alias{twopcor}
\alias{print.twocor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence intervals for two-sided tests on correlation coefficients.
}
\description{
The \code{twopcor} function tests whether the difference between two Pearson correlations is 0. The \code{twocor} function performs the same test on a robust correlation coefficient (percentage bend correlation or Winsorized correlation).
}
\usage{
twopcor(x1, y1, x2, y2, nboot = 599, ...)
twocor(x1, y1, x2, y2, corfun = "pbcor", nboot = 599, tr = 0.2, beta = 0.2, ...)
}
\arguments{
  \item{x1}{
a numeric vector.
}
  \item{y1}{
a numeric vector.
}
  \item{x2}{
a numeric vector.
}
  \item{y2}{
a numeric vector.
}
 \item{nboot}{
number of bootstrap samples.
}
  \item{corfun}{
Either \code{"pbcor"} for percentage based correlation or \code{"wincor"} for Winsorized correlation.
}
 \item{tr}{
amount of Winsorization.
 }
 \item{beta}{
bending constant.
}
\item{...}{
currently ignored.
}
}
\details{
It is tested whether the first correlation coefficient (based on \code{x1} and \code{y1}) equals to the second correlation coefficient (based on \code{x2} and \code{y2}). Both approaches return percentile bootstrap CIs.
}
\value{
\code{twopcor} and \code{twocor} return an object of class \code{"twocor"} containing:

\item{r1}{robust correlation coefficient}
\item{r2}{value of the test statistic}
\item{ci}{confidence interval}
\item{p.value}{p-value}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{pbcor}}, \code{\link{wincor}}
}
\examples{
ct1 <- subset(hangover, subset = (group == "control" & time == 1))$symptoms
ct2 <- subset(hangover, subset = (group == "control" & time == 2))$symptoms
at1 <- subset(hangover, subset = (group == "alcoholic" & time == 1))$symptoms
at2 <- subset(hangover, subset = (group == "alcoholic" & time == 2))$symptoms

set.seed(111)
twopcor(ct1, ct2, at1, at2)
set.seed(123)
twocor(ct1, ct2, at1, at2, corfun = "pbcor", beta = 0.15)
set.seed(224)
twocor(ct1, ct2, at1, at2, corfun = "wincor", tr = 0.15, nboot = 50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

