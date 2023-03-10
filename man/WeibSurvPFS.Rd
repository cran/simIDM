% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalFunctions.R
\name{WeibSurvPFS}
\alias{WeibSurvPFS}
\title{PFS Survival Function from Weibull Transition Hazards}
\usage{
WeibSurvPFS(t, h01, h02, p01, p02)
}
\arguments{
\item{t}{(\code{numeric})\cr  study time-points.}

\item{h01}{(positive \code{number})\cr transition hazard for 0 to 1 transition.}

\item{h02}{(positive \code{number})\cr transition hazard for 0 to 2 transition.}

\item{p01}{(positive \code{number})\cr rate parameter of Weibull distribution for \code{h01}.}

\item{p02}{(positive \code{number})\cr rate parameter of Weibull distribution for \code{h02}.}
}
\value{
This returns the value of PFS survival function at time t.
}
\description{
PFS Survival Function from Weibull Transition Hazards
}
\examples{
WeibSurvPFS(c(1:5), 0.2, 0.5, 1.2, 0.9)
}
