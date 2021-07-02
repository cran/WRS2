\name{trimcibt}
\alias{trimcibt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap-t method for one-sample test}
\description{
Compute a 1-alpha confidence interval for the trimmed mean using a bootstrap percentile t method.
}
\usage{
trimcibt(x, nv = 0, tr = 0.2, alpha = 0.05, nboot = 200, ...)
}

\arguments{
  \item{x}{
    a numeric vector.
  }
  \item{nv}{
    value under H0.
  }
  \item{tr}{
    trim level for the mean.
  }
  \item{alpha}{
    alpha level.
  }
  \item{nboot}{
    number of bootstrap samples.
  }
  \item{...}{
    currently ignored.
  }
}

\value{
  Returns an object of class \code{"trimcibt"} containing:

  \item{ci}{95\% confidence interval}
  \item{estimate}{trimmed mean}
  \item{p.value}{p-value}
  \item{test.stat}{t-statistic}
  \item{tr}{trimming level}
  \item{n}{number of effective observations}
}
\references{
  Wilcox, R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.). Elsevier.
}
\seealso{
  \code{\link{onesampb}}
}
\examples{
  set.seed(123)
  x <- rnorm(30)
  trimcibt(x, nboot = 100)    ## H0: Psi = 0
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

