\name{tsplit}
\alias{tsplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A robust two-way mixed ANOVA (split-plot) for trimmed means.
}
\description{
This function computes a two-way between-within subjects ANOVA on the trimmed means. It allows for one between subjects-variable and for
one within-subjects variable.
}
\usage{
tsplit(formula, random, data, tr = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
 \item{random}{
within subject specification (random effects style)
}
  \item{data}{
an optional data frame for the input data.
}
  \item{tr}{
trim level for the mean.
}
}
\value{
Returns an object of class \code{t2way} containing:

\item{Qa}{first main effect}
\item{A.p.value}{p-value first main effect}
\item{Qb}{second main effect}
\item{B.p.value}{p-value second main effect}
\item{Qab}{interaction effect}
\item{AB.p.value}{p-value interaction effect}
\item{call}{function call}
\item{varnames}{variable names}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{t2way}} 
}
\examples{
## data need to be on long format
pictureLong <- reshape(picture, direction = "long", varying = list(3:4), idvar = "case", 
timevar = c("pictype"), times = c("couple", "alone"))
colnames(pictureLong)[4] <- "friend_requests"

tsplit(friend_requests ~ relationship_status*pictype, random = ~1|case/pictype, data = pictureLong)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

