\name{t3way}
\alias{t3way}
\alias{print.t3way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A three-way ANOVA for trimmed means.
}
\description{
This function computes a three-way ANOVA for trimmed means with all interactions effects. 
}
\usage{
t3way(formula, data, tr = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
  \item{tr}{
trim level for the mean.
}
}
\value{
Returns an object of class \code{t3way} containing:

\item{Qa}{first main effect}
\item{A.p.value}{p-value first main effect}
\item{Qb}{second main effect}
\item{B.p.value}{p-value second main effect}
\item{Qc}{third main effect}
\item{C.p.value}{p-value third main effect}
\item{Qab}{first two-way interaction effect}
\item{AB.p.value}{p-value first two-way interaction effect}
\item{Qac}{second two-way interaction effect}
\item{AC.p.value}{p-value second two-way interaction effect}
\item{Qbc}{third two-way interaction effect}
\item{BC.p.value}{p-value third two-way interaction effect}
\item{Qabc}{three-way interaction effect}
\item{ABC.p.value}{p-value three-way interaction effect}
\item{call}{function call}
\item{varnames}{variable names}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
code{\link{t1way}}, code{\link{t2way}}
}
\examples{
t3way(aggressive ~ degree*gender*type, data = movie)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

