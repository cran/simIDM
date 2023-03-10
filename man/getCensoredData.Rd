% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eventTracking.R
\name{getCensoredData}
\alias{getCensoredData}
\title{Helper function for \code{censoringByNumberEvents}}
\usage{
getCensoredData(time, event, data)
}
\arguments{
\item{time}{(\code{numeric}) \cr event times.}

\item{event}{(\code{numeric})\cr event indicator.}

\item{data}{(\code{data.frame})\cr data frame including patient id \code{id}, recruiting time \code{recruitTime}
and individual censoring time \code{censTimeInd}.}
}
\value{
This function returns a data frame with columns:
event time, censoring indicator, event indicator and event time
in calendar time.
}
\description{
Helper function for \code{censoringByNumberEvents}
}
\examples{
transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

simStudy <- getOneClinicalTrial(
  nPat = c(20, 20), transitionByArm = list(transition1, transition2),
  dropout = list(rate = 0.3, time = 10),
  accrual = list(param = "time", value = 7)
)
simStudyWide <- getDatasetWideFormat(simStudy)
simStudyWide$censTimeInd <- 5 - simStudyWide$recruitTime
NotRecruited <- simStudyWide$id[simStudyWide$censTimeInd < 0]
censoredData <- simStudyWide[!(simStudyWide$id \%in\% NotRecruited), ]
getCensoredData(time = censoredData$OStime, event = censoredData$OSevent, data = censoredData)
}
