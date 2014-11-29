\name{rmanova}
\alias{rmanova}
\alias{rmmcp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A heteroscedastic one-way repeated measures ANOVA for trimmed means.
}
\description{
The \code{rmanova} function computes a one-way repeated measures ANOVA for the trimmed means. Homoscedasticity assumption not required.  
Corresponding post hoc tests can be performed using \code{rmmcp}.
}
\usage{
rmanova(y, groups, blocks, tr = 0.2)
rmmcp(y, groups, blocks, tr = 0.2)
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
}
\value{
\code{rmanova} an object of class \code{"t1way"} containing:

\item{test}{value of the test statistic}
\item{df1}{degrees of freedom}
\item{df2}{degrees of freedom}
\item{p.value}{p-value}
\item{call}{function call}

\code{rmmcp} returns an object of class \code{"mcp1"} containing:

\item{comp}{inference for all pairwise comparisons}
\item{fnames}{names of the factor levels}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{med1way}},\code{\link{t1way}}
}
\examples{
head(WineTasting)
rmanova(WineTasting$Taste, WineTasting$Wine, WineTasting$Taster)

## post hoc
rmmcp(WineTasting$Taste, WineTasting$Wine, WineTasting$Taster)

head(bush)
require(reshape)
bushLong <- melt(bush, id.var = "participant", variable_name = "food")
rmanova(bushLong$value, bushLong$food, bushLong$participant)

## post hoc
rmmcp(bushLong$value, bushLong$food, bushLong$participant)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

