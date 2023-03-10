% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSimulatedData.R
\name{addStaggeredEntry}
\alias{addStaggeredEntry}
\title{Staggered Study Entry}
\usage{
addStaggeredEntry(simData, N, accrualParam, accrualValue)
}
\arguments{
\item{simData}{(\code{data.frame})\cr simulated data frame containing entry and exit times
at individual study time scale. See \code{\link[=getSimulatedData]{getSimulatedData()}} for details.}

\item{N}{(\code{int})\cr number of patients.}

\item{accrualParam}{(\code{string})\cr possible values are 'time' or 'intensity'.}

\item{accrualValue}{(\code{number})\cr specifies the accrual intensity. For \code{accrualParam} equal time,
it describes the length of the accrual period. For \code{accrualParam} equal intensity, it describes
the number of patients recruited per time unit.  If \code{accrualValue} is equal to 0,
all patients start at calendar time 0
in the initial state.}
}
\value{
This returns a data set containing a single simulated study containing accrual times,
i.e. staggered study entry.
This is a helper function of \code{\link[=getSimulatedData]{getSimulatedData()}}.
}
\description{
This function adds staggered study entry times to a simulated data set with illness-death model structure.
}
\examples{
simData <- data.frame(
  id = c(1, 1, 2, 3), from = c(0, 1, 0, 0), to = c(1, 2, "cens", 2),
  entry = c(0, 3, 0, 0),
  exit = c(3, 5.3, 5.6, 7.2), censTime = c(6.8, 6.8, 5.6, 9.4)
)
addStaggeredEntry(simData, 3, accrualParam = "time", accrualValue = 5)
}
