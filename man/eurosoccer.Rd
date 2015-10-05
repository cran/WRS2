\name{eurosoccer}
\alias{eurosoccer}
\docType{data}
\title{European Soccer Leagues}
\description{Contains various team stats from five European soccer leagues (2008/09 season).
}
\usage{eurosoccer}
\format{
  A data frame with 11 variables and 96 teams:
  \describe{
     \item{\code{League}}{Country}                  
     \item{\code{Team}}{Team} 
     \item{\code{Games}}{Number of games} 
     \item{\code{Won}}{Games won} 
     \item{\code{Tied}}{Games tied} 
     \item{\code{Lost}}{Games lost} 
     \item{\code{GoalsScored}}{Goals scored} 
     \item{\code{GoalsConceded}}{Goals conceded} 
     \item{\code{GoalDifference}}{Goal difference} 
     \item{\code{Points}}{Final amount of points} 
     \item{\code{GoalsGame}}{Goal scored per game}
  }
}

\examples{
head(eurosoccer)
}
\keyword{datasets}
