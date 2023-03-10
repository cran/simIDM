% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getClinicalTrials.R
\name{getOneClinicalTrial}
\alias{getOneClinicalTrial}
\title{Simulation of a Single Oncology Clinical Trial}
\usage{
getOneClinicalTrial(
  nPat,
  transitionByArm,
  dropout = list(rate = 0, time = 12),
  accrual = list(param = "time", value = 0)
)
}
\arguments{
\item{nPat}{(\code{integer})\cr numbers of patients per treatment arm.}

\item{transitionByArm}{(\code{list}) \cr transition parameters for each treatment group.
See \code{\link[=exponential_transition]{exponential_transition()}}, \code{\link[=piecewise_exponential]{piecewise_exponential()}} and \code{\link[=weibull_transition]{weibull_transition()}} for details.}

\item{dropout}{dropout (\code{list})\cr specifies drop-out probability. See \code{\link[=getSimulatedData]{getSimulatedData()}} for details.
Can be specified either as one list that should be applied to all treatment groups or a separate list
for each treatment group.}

\item{accrual}{accrual (\code{list})\cr specifies accrual intensity. See \code{\link[=addStaggeredEntry]{addStaggeredEntry()}} for details.
Can be specified either as one list that should be applied to all treatment groups or a separate list
for each treatment group.}
}
\value{
This returns a data frame with one simulated clinical trial and multiple treatment arms.
See \code{\link[=getSimulatedData]{getSimulatedData()}} for the explanation of the columns. The column \code{trt} contains the treatment indicator.
This is a helper function of \code{\link[=getClinicalTrials]{getClinicalTrials()}}.
}
\description{
This function creates a data set with a single simulated oncology clinical trial with one row per transition
based on an illness-death model. Studies with an arbitrary number of treatment arms are possible.
}
\examples{
transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
transition3 <- exponential_transition(h01 = 1.1, h02 = 1, h12 = 1.5)
getOneClinicalTrial(
  nPat = c(30, 20, 30), transitionByArm = list(transition1, transition2, transition3),
  dropout = list(rate = 0, time = 12),
  accrual = list(param = "time", value = 0)
)
}
