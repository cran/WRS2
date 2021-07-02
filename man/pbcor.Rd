\name{pbcor}
\alias{pbcor}
\alias{wincor}
\alias{winall}
\alias{pball}
\alias{print.pbcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust correlation coefficients.
}
\description{
The \code{pbcor} function computes the percentage bend correlation coefficient, \code{wincor} the Winsorized correlation,
\code{pball} the percentage bend correlation matrix, \code{winall} the Winsorized correlation matrix.
}
\usage{
pbcor(x, y = NULL, beta = 0.2, ci = FALSE, nboot = 500, alpha = 0.05, ...)
pball(x, beta = 0.2, ...)
wincor(x, y = NULL, tr = 0.2, ci = FALSE, nboot = 500, alpha = 0.05, ...)
winall(x, tr = 0.2, ...)
}
\arguments{
  \item{x}{
a numeric vector, a matrix or a data frame.
}
 \item{y}{
a second numeric vector (for correlation functions).
}
  \item{beta}{
bending constant.
}
 \item{tr}{
amount of Winsorization.
 }
\item{ci}{
whether boostrap CI should be computed or not.
 }
  \item{nboot}{
number of bootstrap samples for CI computation.
 }
  \item{alpha}{
alpha level for CI computation.
 }
 \item{...}{
currently ignored.
}
}
\details{
It tested is whether the correlation coefficient equals 0 (null hypothesis) or not. Missing values are deleted pairwise. The tests are sensitive to heteroscedasticity. The test statistic H in \code{pball} tests the hypothesis that all correlations are equal to zero.
}
\value{
\code{pbcor} and \code{wincor} return an object of class \code{"pbcor"} containing:

\item{cor}{robust correlation coefficient}
\item{test}{value of the test statistic}
\item{p.value}{p-value}
\item{n}{number of effective observations}
\item{cor_ci}{bootstrap confidence interval}
\item{call}{function call}

\code{pball} and \code{winall} return an object of class \code{"pball"} containing:

\item{pbcorm}{robust correlation matrix}
\item{p.values}{p-values}
\item{H}{H-statistic}
\item{H.p.value}{p-value H-statistic}
\item{cov}{variance-covariance matrix}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{twocor}}
}
\examples{
x1 <- subset(hangover, subset = (group == "control" & time == 1))$symptoms
x2 <- subset(hangover, subset = (group == "control" & time == 2))$symptoms

pbcor(x1, x2)
pbcor(x1, x2, beta = 0.1, ci = TRUE)

wincor(x1, x2)
wincor(x1, x2, tr = 0.1, ci = TRUE)

require(reshape)
hanglong <- subset(hangover, subset = group == "control")
hangwide <- cast(hanglong, id ~ time, value = "symptoms")[,-1]

pball(hangwide)
winall(hangwide)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

