\name{trimse}
\alias{trimse}
\alias{winmean}
\alias{winse}
\alias{msmedse}
\alias{mest}
\alias{mestse}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust location measures and their standard errors (se).
}
\description{
The following functions for estimating robust location measures and their standard errors are provided: \code{winmean} 
for the Winsorized mean, \code{winse} for its se, \code{trimse} for the trimmend mean se, \code{msmedse} for the median se, 
\code{mest} for the M-estimator with se in \code{mestse}. 
}
\usage{
winmean(x, tr = 0.2, na.rm = FALSE)
winse(x, tr = 0.2)
trimse(x, tr = 0.2, na.rm = FALSE)
msmedse(x, sewarn = TRUE)
mest(x, bend = 1.28, na.rm = FALSE)
mestse(x, bend = 1.28)
}
\arguments{
  \item{x}{
a numeric vector containing the values whose measure is to be computed.
}
  \item{tr}{
trim lor Winsorizing level. 
}
\item{na.rm}{
a logical value indicating whether NA values should be stripped before the computation proceeds.
}
\item{sewarn}{
a logical value indicating whether warnings for ties should be printed.
}
\item{bend}{
bending constant for M-estimator.
}
}

\details{
The standard error for the median is computed according to McKean and Shrader (1984).
}

\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.

McKean, J. W., & Schrader, R. M. (1984). A comparison of methods for studentizing the sample median. 
Communications in Statistics - Simulation and Computation, 13, 751-773. 

Dana, E. (1990). Salience of the self and salience of standards: Attempts to match self to standard. 
Unpublished PhD thesis, Department of Psychology, University of Southern California. 
}

\examples{
## Self-awareness data (Dana, 1990): Time persons could keep a portion of an 
## apparatus in contact with a specified range.
self <- c(77, 87, 88, 114, 151, 210, 219, 246, 253, 262, 296, 299, 306, 376, 
          428, 515, 666, 1310, 2611)
mean(self, 0.1)     ## .10 trimmed mean 
trimse(self, 0.1)   ## se trimmed mean
winmean(self, 0.1)  ## Winsorized mean (.10 Winsorizing amount)
winse(self, 0.1)    ## se Winsorized mean
median(self)        ## median
msmedse(self)       ## se median
mest(self)          ## Huber M-estimator
mestse(self) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

