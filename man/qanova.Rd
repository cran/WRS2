\name{Qanova}
\alias{Qanova}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantile ANOVA
}
\description{
One-way ANOVA based on quantiles. Only known method to work well when tied values are likely to occur.
}
\usage{
Qanova(formula, data, q = 0.5, nboot = 600, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
  \item{q}{
quantile (or vector of quantiles) to be used.
}
\item{nboot}{
number of bootstrap samples
}
\item{...}{
currently ignored.
}
}
\details{
Test global hypothesis that J independent groups have equal quantiles (default: median) using the Harrell-Davis estimator. Performs well when there are tied values.
}
\value{
\code{Qanova} an object of class \code{"qanova"} containing:

\item{psihat}{value of the test statistics}
\item{contrasts}{contrasts}
\item{p.value}{p-value}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{qcomhd}}, \code{\link{t1way}}, \code{\link{discANOVA}}
}
\examples{
## median comparison
set.seed(123)
fitmed <- Qanova(libido ~ dose, viagra, nboot = 200)
fitmed

## 1st, 3rd quartile comparison
set.seed(123)
fitquart <- Qanova(libido ~ dose, viagra, q = c(0.25, 0.75), nboot = 200)
fitquart



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

