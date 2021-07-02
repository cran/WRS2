\name{wmcpAKP}
\alias{wmcpAKP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effect size for dependent samples ANOVA
}
\description{
Compute an AKP-type effect size for dependent sample ANOVA
}
\usage{
wmcpAKP(x, tr = 0.2, nboot = 200, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data frame in wide format (no missing values allowed)
}
  \item{tr}{
trim level for the means.
}
  \item{nboot}{
number of bootstrap samples.
}
\item{...}{
currently ignored.
}
}
\details{
The computation is based on a modification of the Algina-Keselman-Penfield effect size for J dependent samples.
}


\references{
Algina, J., Keselman, H.J., & Penfield, R.D. (2005). An alternative to Cohen's standardized mean difference effect size: A robust parameter and confidence interval in the two independent groups case. Psychological Methods, 10, 317-328.

}
\seealso{
\code{\link{rmanova}}, \code{\link{yuend}}
}

\examples{
\dontrun{
require(reshape)
WineLong <- cast(WineTasting, Taster ~ Wine, value = "Taste")[,-1]
set.seed(123)
effsize <- wmcpAKP(WineLong, nboot = 20)
effsize
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

