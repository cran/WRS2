\name{yuend}
\alias{yuend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Paired samples $t$-test on trimmed means.
}
\description{
The function \code{yuen} performs Yuen's test for trimmed means, \code{yuenbt} is a bootstrap version of it. The \code{pb2gen} function performs a t-test based on various robust estimators. 
}
\usage{
yuend(x, y, tr = 0.2)
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

}
\details{
The test statistic is a paired samples generalization of Yuen's independent samples t-test on trimmed means.
}

\value{
Returns an object of class \code{"yuen"} containing:

\item{test}{value of the test statistic (t-statistic)}
\item{p.value}{p-value}
\item{conf.int}{confidence interval}
\item{df}{degress of freedom}
\item{diff}{trimmed mean difference}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{yuen}}
}

\examples{
## Cholesterol data from Wilcox (2012, p. 197)
before <- c(190, 210, 300,240, 280, 170, 280, 250, 240, 220)
after <- c(210, 210, 340, 190, 260, 180, 200, 220, 230, 200)
yuend(before, after)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

