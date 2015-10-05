\name{runmean}
\alias{runmean}
\alias{rungen}
\alias{runmbo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Running Interval Smoother
}
\description{
The \code{runmean} implements a running interval smoother on the trimmed mean, \code{rungen} uses general M-estimators, \code{runmbo} performs 
interval smoothing on M-estimators with bagging. 
}
\usage{
runmean(x, y, fr = 1, tr = 0.2)
rungen(x, y, fr = 1, est = "mom")
runmbo(x, y, fr = 1, est = "mom", nboot = 40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector of data values (predictor)
}
  \item{y}{
a numeric vector of data values (response)
}
\item{fr}{
smoothing factor (see details)
}
  \item{tr}{
trim level for the mean
}
  \item{est}{
type of M-estimator (\code{"mom"}, \code{"onestep"}, or \code{"median"}) 
}
  \item{nboot}{
number of bootstrap samples
}
}

\details{
The larger the smoothing factor, the stronger the smoothing. Often the choice \code{fr = 1} gives good results; the general strategy is 
to find the smallest constant so that the plot looks reasonably smooth. 
}
\value{
Returns the fitted values.
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.
}
\seealso{
\code{\link{ancova}}
}
\examples{
## trimmed mean smoother
fitmean <- runmean(Pygmalion$Pretest, Pygmalion$Posttest)
## MOM smoother
fitmest <- rungen(Pygmalion$Pretest, Pygmalion$Posttest)
## median smoother
fitmed <- rungen(Pygmalion$Pretest, Pygmalion$Posttest, est = "median")
## bagged onestep smoother
fitbag <- runmbo(Pygmalion$Pretest, Pygmalion$Posttest, est = "onestep")

## plot smoothers
plot(Pygmalion$Pretest, Pygmalion$Posttest, col = "gray", xlab = "Pretest", ylab = "Posttest",
     main = "Pygmalion Smoothing")
orderx <- order(Pygmalion$Pretest)
lines(Pygmalion$Pretest[orderx], fitmean[orderx], lwd = 2)
lines(Pygmalion$Pretest[orderx], fitmest[orderx], lwd = 2, col = 2)
lines(Pygmalion$Pretest[orderx], fitmed[orderx], lwd = 2, col = 3)
lines(Pygmalion$Pretest[orderx], fitbag[orderx], lwd = 2, col = 4)
legend("topleft", legend = c("Trimmed Mean", "MOM", "Median", "Bagged Onestep"), col = 1:4, lty = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

