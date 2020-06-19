\name{yuend}
\alias{yuend}
\alias{Dqcomhd}
\alias{dep.effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Paired samples robust t-tests.
}
\description{
The function \code{yuend} performs Yuen's test on trimmed means for dependent samples. \code{Dqcomhd} compares the quantiles of the marginal distributions associated with two dependent groups via hd estimator. Tied values are allowed.
\code{dep.effect} computes various effect sizes and confidence intervals for two dependent samples (see Details). 
}
\usage{
yuend(x, y, tr = 0.2)
Dqcomhd(x, y, q = c(1:9)/10, nboot = 1000, na.rm = TRUE)
dep.effect(x, y, tr = 0.2, nboot = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an numeric vector of data values (e.g. for time 1).
}
 \item{y}{
an numeric vector of data values (e.g. for time 2).
}
  \item{tr}{
trim level for the means.
}
  \item{q}{
quantiles to be compared.
}
  \item{nboot}{
number of bootstrap samples.
}
  \item{na.rm}{
whether missing values should be removed. 
}

}
\details{
The test statistic is a paired samples generalization of Yuen's independent samples t-test on trimmed means.

\code{dep.effect} computes the following effect sizes: 

AKP: robust standardized difference similar to  Cohen's d

QS:  Quantile shift based on the median of the distribution of difference scores, 

QStr: Quantile shift based on the trimmed mean of the distribution of X-Y 

SIGN:  P(X<Y), probability that for a random pair, the first is less than the second.
}

\value{
\code{yuend} returns an object of class \code{"yuen"} containing:

\item{test}{value of the test statistic (t-statistic)}
\item{p.value}{p-value}
\item{conf.int}{confidence interval}
\item{df}{degress of freedom}
\item{diff}{trimmed mean difference}
\item{call}{function call}

\code{Dqcomhd} returns an object of class \code{"robtab"} containing:

\item{partable}{parameter table}

\code{dep.effect} returns a matrix with the null value of the effect size, the estimated effect size, small/medium/large conventions, and lower/upper CI bounds.  

}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{yuen}}, \code{\link{qcomhd}} 
}

\examples{
## Cholesterol data from Wilcox (2012, p. 197)
before <- c(190, 210, 300,240, 280, 170, 280, 250, 240, 220)
after <- c(210, 210, 340, 190, 260, 180, 200, 220, 230, 200)
yuend(before, after)

set.seed(123)
Dqcomhd(before, after, nboot = 200, q = c(0.25, 0.5, 0.75))

set.seed(123)
dep.effect(before, after)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

