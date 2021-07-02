\name{ZYmediate}
\alias{ZYmediate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust mediation test
}
\description{
Performs a robust mediation test according to the method proposed by Zu & Yuan (2010).
}
\usage{
ZYmediate(x, y, med, nboot = 2000, alpha = 0.05, kappa = 0.05, ...)
}
\arguments{
  \item{x}{
vector with predictor values.
}
  \item{y}{
vector with response values.
}
 \item{med}{
vector with mediator values.
}
  \item{nboot}{
number of bootstrap samples.
}
\item{alpha}{
alpha level.
}
\item{kappa}{
the percent of cases to be controlled when robust method is used.
}
\item{...}{
currently ignored.
}
}

\value{
\code{ZYmediate} returns an object of class \code{"robmed"} containing:

\item{a.est}{effect of predictor on mediator)}
\item{b.est}{effect of mediator on response (in multiple regression model which includes the predictor as well)}
\item{ab.est}{indirect effect (mediation effect)}
\item{CI.ab}{confidence interval mediation effect}
\item{p.value}{p-value mediation test}
\item{call}{function call}

}
\references{
Zu J., Yuan, K.-H. (2010). Local influence and robust procedures for mediation analysis. Multivariate Behavioral Research, 45, 1-44.
}
\examples{
## Leerkes data:
## Y: Efficacy
## X: MatCare
## M: Esteem

## fitting robust mediator regressions
require(MASS)
summary(rlm(Efficacy ~ MatCare, data = Leerkes))
summary(rlm(Esteem ~ MatCare, data = Leerkes))
summary(rlm(Efficacy ~ MatCare + Esteem, data = Leerkes))

## robust testing of mediating effect (indirect effect)
with(Leerkes, ZYmediate(MatCare, Efficacy, Esteem))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

