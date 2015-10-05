\name{swimming}
\alias{swimming}
\docType{data}
\title{Optimistic and Pessimistic Swimmers}
\description{At a swimming team practice, all participants were asked to swim their best event as far as possible, but in each case the time that was reported was falsified to indicate poorer than expected performance (i.e., each swimmer was disappointed). 30 min later, they did the same performance. The authors predicted that on the second trial more pessimistic swimmers would do worse than on their first trial,
whereas optimistic swimmers would do better. The response is ratio = Time1/Time2 (> 1 means that a swimmer did better in trial 2).
}
\usage{swimming}
\format{
  A data frame with 4 variables and 58 observations:
  \describe{
    \item{\code{Optim}}{Optimists and pessimists}
    \item{\code{Sex}}{Gender of the swimmer}
    \item{\code{Event}}{Swimming event: freestyle, breaststroke, backstroke}
    \item{\code{Ratio}}{Ratio between the swimming times}
  }
}

\references{
Seligman, M. E. P., Nolen-Hoeksema, S., Thornton, N., & Thornton, C. M. (1990). Explanatory style as a mechanism of disappointing athletic 
performance. Psychological Science, 1, 143-146. 
}
\examples{
summary(swimming)
}
\keyword{datasets}
