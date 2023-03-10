% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getClinicalTrials.R
\name{getDatasetWideFormat}
\alias{getDatasetWideFormat}
\title{Conversion of a Data Set from One Row per Transition to One Row per Patient}
\usage{
getDatasetWideFormat(data)
}
\arguments{
\item{data}{(\code{data.frame})\cr data frame containing entry and exit times of an illness-death model.
See \code{\link[=getSimulatedData]{getSimulatedData()}} for details.}
}
\value{
This function returns a data set with one row per patient and endpoints PFS and OS.
}
\description{
Conversion of a Data Set from One Row per Transition to One Row per Patient
}
\details{
The output data set contains the following columns:
\itemize{
\item id (\code{integer}): patient id.
\item trt \code{integer}): treatment id.
\item PFStime (\code{numeric}): event time of PFS event.
\item CensoredPFS (\code{logical}): censoring indicator for PFS event.
\item PFSevent (\code{logical}): event indicator for PFS event.
\item OStime (\code{numeric}): event time of OS event.
\item CensoredOS (\code{logical}): censoring indicator for OS event.
\item OSevent (\code{logical}): event indicator for OS event.
\item recruitTime (\code{numeric}): time of recruitment.
\item OStimeCal (\code{numeric}): OS event time at calendar time scale.
\item PFStimeCal (\code{numeric}): PFS event time at calendar time scale.
}
}
\examples{
transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
transition3 <- exponential_transition(h01 = 1.1, h02 = 1, h12 = 1.5)
simData <- getOneClinicalTrial(
  nPat = c(30, 20, 30), transitionByArm = list(transition1, transition2, transition3),
  dropout = list(rate = 0, time = 12),
  accrual = list(param = "time", value = 0)
)
getDatasetWideFormat(simData)
}
