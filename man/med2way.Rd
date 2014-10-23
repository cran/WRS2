\name{med2way}
\alias{med2way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A two-way ANOVA for medians.
}
\description{
This function computes a two-way ANOVA medians with interactions effects. 
}
\usage{
med2way(formula, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
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
\code{\link{t2way}}, \code{\link{med1way}} 
}
\examples{
med2way(attractiveness ~ gender*alcohol, data = goggles)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

