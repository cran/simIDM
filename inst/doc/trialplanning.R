## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE,  fig.cap = " Figure 1 - Multistate model with indermediate state progession and absorbing state death",fig.height = 5, fig.width = 8, out.width = "60%"----
library(prodlim)
plotIllnessDeathModel(
  style = 1, box1.label = "0: initial state", box2.label = "1: progression",
  box3.label = "2: death", arrowLabelSymbol = "h"
)

## -----------------------------------------------------------------------------
library(simIDM)
library(survival)
transitionTrt <- exponential_transition(h01 = 0.2, h02 = 0.28, h12 = 0.4)
transitionPl <- exponential_transition(h01 = 0.4, h02 = 0.3, h12 = 0.5)

transitionList <- list(transitionPl, transitionTrt)

## -----------------------------------------------------------------------------
timepoints <- c(0, 0.1, 0.2, 0.3, 0.7, 1, 5)
ExpSurvOS(timepoints, h01 = 0.2, h02 = 0.4, h12 = 0.1)
WeibSurvOS(timepoints, h01 = 0.2, h02 = 0.5, h12 = 2.1, p01 = 1.2, p02 = 0.9, p12 = 1)
PWCsurvOS(timepoints,
  h01 = c(0.3, 0.5), h02 = c(0.5, 0.8), h12 = c(0.7, 1),
  pw01 = c(0, 4), pw02 = c(0, 8), pw12 = c(0, 3)
)

## -----------------------------------------------------------------------------
transitionListNULL <- list(transitionPl, transitionPl)
nRep <- 100
SimNULL <- getClinicalTrials(
  nRep = nRep, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionListNULL,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 100)
)

## -----------------------------------------------------------------------------
alphaOS <- 0.04
alphaPFS <- 0.01
CriticalOS <- abs(qnorm(alphaOS / 2))
CriticalPFS <- abs(qnorm(alphaPFS / 2))

## -----------------------------------------------------------------------------
library(coxphw)
avergHR <- function(data) {
  fit <- coxphw(formula = Surv(OStime + 0.001, OSevent) ~ trt, data = data, template = "AHR")
  exp(coef(fit))
}

StudyforAvg <- getClinicalTrials(
  nRep = 100, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionList,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 25)
)

averageHR <- sapply(StudyforAvg, avergHR)
MeanAverageHR <- mean(averageHR)

## -----------------------------------------------------------------------------
LogRankTest <- function(data, endpoint, critical) {
  time <- if (endpoint == "OS") {
    data$OStime
  } else if (endpoint == "PFS") {
    data$PFStime
  }
  event <- if (endpoint == "OS") {
    data$OSevent
  } else if (endpoint == "PFS") {
    data$PFSevent
  }
  LogRank <- survdiff(Surv(time, event) ~ trt, data)
  Passed <- sqrt(LogRank$chisq) > critical
  return(Passed)
}

## -----------------------------------------------------------------------------
# Step 1: cut simulated data at time-point of OS/PFS analysis
studyAtPFSAna <- lapply(SimNULL, censoringByNumberEvents,
  eventNum = 329, typeEvent = "PFS"
)

studyAtOSAna <- lapply(SimNULL, censoringByNumberEvents,
  eventNum = 660, typeEvent = "OS"
)

# Step 2: get results of log-rank test for both endpoints and all studies

TestPassedPFS <- lapply(studyAtPFSAna, LogRankTest, endpoint = "PFS", CriticalPFS)
TestPassedOS <- lapply(studyAtOSAna, LogRankTest, endpoint = "OS", CriticalOS)

## -----------------------------------------------------------------------------
# empirical type-I error of PFS
empAlphaPFS <- 100 * (sum(unlist(TestPassedPFS)) / nRep)

# empirical type-I error of OS
empAlphaOS <- 100 * (sum(unlist(TestPassedOS)) / nRep)

TestBoth <- (unlist(TestPassedPFS) + unlist(TestPassedOS) >= 1)
# empirical global significance level
empGlobalAlpha <- 100 * (sum(TestBoth) / nRep)

## ----eval=FALSE---------------------------------------------------------------
#  while (empGlobalAlpha < 4.9) {
#    CriticalOS <- CriticalOS - 0.01
#    CriticalPFS <- CriticalPFS - 0.01
#  
#    TestPassedPFS <- lapply(studyAtPFSAna, LogRankTest, endpoint = "PFS", CriticalPFS)
#    TestPassedOS <- lapply(studyAtOSAna, LogRankTest, endpoint = "OS", CriticalOS)
#    TestBoth <- (unlist(TestPassedPFS) + unlist(TestPassedOS) >= 1)
#  
#    # empirical global significance level
#    empGlobalAlpha <- 100 * (sum(TestBoth) / nRep)
#  }

## -----------------------------------------------------------------------------

SimH1 <- getClinicalTrials(
  nRep = nRep, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionList,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 100)
)

## -----------------------------------------------------------------------------
# Step 1: cut simulated data at time-point of OS/PFS analysis
studyAtPFSAnaH1 <- lapply(SimH1, censoringByNumberEvents,
  eventNum = 329, typeEvent = "PFS"
)

studyAtOSAnaH1 <- lapply(SimH1, censoringByNumberEvents,
  eventNum = 660, typeEvent = "OS"
)

# Step 2: get results of log-rank test for both endpoints and all studies.

logrankPFSH1 <- lapply(studyAtPFSAnaH1, LogRankTest,
  endpoint = "PFS", CriticalPFS
)
logrankOSH1 <- lapply(studyAtOSAnaH1, LogRankTest,
  endpoint = "OS", CriticalOS
)

TestPassedPFSH1 <- lapply(logrankPFSH1, `[[`, 1)
TestPassedOSH1 <- lapply(logrankOSH1, `[[`, 1)
# Step 3: count significant log-rank tests.
# empirical power PFS
empPowerPFS <- 100 * (sum(unlist(TestPassedPFSH1)) / nRep)
empPowerPFS
# empirical power OS
empPowerOS <- 100 * (sum(unlist(TestPassedOSH1)) / nRep)
empPowerOS
# joint power
TestBothH1 <- (unlist(TestPassedPFSH1) + unlist(TestPassedOSH1) == 2)

jointPower <- 100 * (sum(TestBothH1) / nRep)
jointPower

## -----------------------------------------------------------------------------
# median time
TimePointsPFS <- lapply(SimH1, getTimePoint,
  eventNum = 329, typeEvent = "PFS",
  byArm = FALSE
)
median_timePFS <- median(unlist(TimePointsPFS))

TimePointsOS <- lapply(SimH1, getTimePoint,
  eventNum = 684, typeEvent = "OS",
  byArm = FALSE
)
median_timeOS <- median(unlist(TimePointsOS))

median_timePFS
median_timeOS

# number of PFS events at time of OS analysis
eventsPFS <- lapply(
  seq_along(TimePointsPFS),
  function(t) {
    return(sum(SimH1[[t]]$OSevent[(SimH1[[t]]$OStime + SimH1[[t]]$recruitTime)
    <= TimePointsPFS[[t]]]))
  }
)
# number of OS events at time of PFS analysis
eventsOS <- lapply(
  seq_along(TimePointsOS),
  function(t) {
    return(sum(SimH1[[t]]$PFSevent[(SimH1[[t]]$PFStime + SimH1[[t]]$recruitTime)
    <= TimePointsOS[[t]]]))
  }
)

nPFS <- mean(unlist(eventsPFS))
nOS <- mean(unlist(eventsOS))

nPFS
nOS

