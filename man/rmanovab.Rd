\name{rmanovab}
\alias{rmanovab}
\alias{pairdepb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A heteroscedastic one-way repeated measures bootstrap ANOVA for trimmed means.
}
\description{
The \code{rmanova} function computes a bootstrap version of the one-way repeated measures ANOVA for the trimmed means. Homoscedasticity assumption not required. Corresponding post hoc tests can be performed using \code{pairdepb}.   
}
\usage{
rmanovab(y, groups, blocks, tr = 0.2, nboot = 599)
pairdepb(y, groups, blocks, tr = 0.2, nboot = 599)
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
}
\value{
\code{rmanovab} an object of class \code{"rmanovab"} containing:

\item{test}{value of the test statistic}
\item{crit}{critical value}
\item{call}{function call}

\code{pairdepb} returns an object of class \code{"mcp2"} containing:

\item{comp}{inference for all pairwise comparisons}
\item{fnames}{names of the factor levels}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{rmanova}}
}
\examples{
head(WineTasting)
rmanovab(WineTasting$Taste, WineTasting$Wine, WineTasting$Taster, nboot = 300)

## post hoc
pairdepb(WineTasting$Taste, WineTasting$Wine, WineTasting$Taster, nboot = 300)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

