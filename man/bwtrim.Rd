\name{bwtrim}
\alias{bwtrim}
\alias{tsplit}
\alias{sppba}
\alias{sppbb}
\alias{sppbi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A robust two-way mixed ANOVA using trimmed means.
}
\description{
The \code{bwtrim} function computes a two-way between-within subjects ANOVA on the trimmed means. It allows for one between subjects variable and one within subjects variable. The functions \code{sppba}, \code{sppbb}, and \code{sppbi} compute the main fixed effect, the main 
within-subjects effect, and the interaction effect only, respectively, using bootstrap. For these 3 functions the user can choose 
an M-estimator for group comparisons. 
}
\usage{
bwtrim(formula, id, data, tr = 0.2)
tsplit(formula, id, data, tr = 0.2)
sppba(formula, id, data, est = "mom", avg = TRUE, nboot = 500, MDIS = FALSE)
sppbb(formula, id, data, est = "mom", nboot = 500)
sppbi(formula, id, data, est = "mom", nboot = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
 \item{id}{
subject ID.
}
  \item{data}{
an optional data frame for the input data.
}
  \item{tr}{
trim level for the mean.
}
 \item{est}{
Estimate to be used for the group comparisons: either \code{"onestep"} for one-step M-estimator of location using Huber's Psi, \code{"mom"} for the modified one-step (MOM) estimator of location based on Huber's Psi, or \code{"median"}. 
}
 \item{avg}{
If \code{TRUE}, the analysis is done by averaging K measures of location for each level of the fixed effect, and then comparing averages by testing the hypothesis that all pairwise differences are equal to zero. If \code{FALSE} the analysis is done by testing whether K equalities are simultaneously true. 
}
 \item{nboot}{
number of bootstrap samples. 
}
 \item{MDIS}{
if \code{TRUE} the depths of the points in the bootstrap cloud are based on Mahalanobis distance, if \code{FALSE} a projection distance is used. 
}
}

\details{
The \code{tsplit} function is doing exactly the same as \code{bwtrim}. It is kept in the package in order to be consistent with 
older versions of the Wilcox (2012) book. 

For \code{sppba},  \code{sppbb}, and \code{sppbi} the analysis is done based on all pairs of difference scores. The null hypothesis is that all such differences have a typical value of zero. In the formula interface it is required to provide full model. 
}

\value{
\code{bwtrim} returns an object of class \code{"t2way"} containing:

\item{Qa}{first main effect}
\item{A.p.value}{p-value first main effect}
\item{Qb}{second main effect}
\item{B.p.value}{p-value second main effect}
\item{Qab}{interaction effect}
\item{AB.p.value}{p-value interaction effect}
\item{call}{function call}
\item{varnames}{variable names}

\code{sppba}, \code{sppbb}, and \code{sppbi} returns an object of class \code{"spp"} containing:

\item{test}{value of the test statistic}
\item{p.value}{p-value}
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

## 2-way within-between subjects ANOVA
bwtrim(friend_requests ~ relationship_status*pictype, id = case, data = pictureLong)

## between groups effect only (MOM estimator)
sppba(friend_requests ~ relationship_status*pictype, case, data = pictureLong)

## within groups effect only (MOM estimator)
sppbb(friend_requests ~ relationship_status*pictype, case, data = pictureLong)

## interaction effect only (MOM estimator)
sppbi(friend_requests ~ relationship_status*pictype, case, data = pictureLong)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

