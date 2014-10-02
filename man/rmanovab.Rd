\name{rmanovab}
\alias{rmanovab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A heteroscedastic one-way repeated measures bootstrap ANOVA for trimmed means.
}
\description{
This function computes a bootstrap version of the one-way repeated measures ANOVA for the trimmed means. Homoscedasticity assumption not required.  
}
\usage{
rmanovab(y, groups, blocks, tr = 0.2, nboot = 599, SEED = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a numeric vector of data values (response).
}
  \item{groups}{
a vector giving the group of the corresponding elements of y.
}
\item{blocks}{
a vector giving the block of the corresponding elements of y.
}
  \item{tr}{
trim level for the mean.
}
 \item{nboot}{
number of bootstrap samples.
}
 \item{SEED}{
whether random number seed should be set.
}
}
\value{
Returns an object of class \code{rmanovab} containing:

\item{test}{value of the test statistic}
\item{crit}{critical value}
\item{call}{function call}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{rmanova}}
}
\examples{
head(WineTasting)
rmanovab(WineTasting$Taste, WineTasting$Wine, WineTasting$Taster)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

