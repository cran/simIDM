% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalFunctions.R
\name{ExpSurvOS}
\alias{ExpSurvOS}
\title{OS Survival Function from Constant Transition Hazards}
\usage{
ExpSurvOS(t, h01, h02, h12)
}
\arguments{
\item{t}{(\code{numeric})\cr  study time-points.}

\item{h01}{(positive \code{number})\cr transition hazard for 0 to 1 transition.}

\item{h02}{(positive \code{number})\cr transition hazard for 0 to 2 transition.}

\item{h12}{(positive \code{number})\cr transition hazard for 1 to 2 transition.}
}
\value{
This returns the value of OS survival function at time t.
}
\description{
OS Survival Function from Constant Transition Hazards
}
\examples{
ExpSurvOS(c(1:5), 0.2, 0.4, 0.1)
}
