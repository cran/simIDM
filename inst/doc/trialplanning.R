## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE,  fig.cap = " Figure 1 - Multistate model with indermediate state progession and absorbing state death", fig.height = 5, fig.width = 8, out.width = "60%"----
library(prodlim)
plotIllnessDeathModel(
  style = 1, box1.label = "0: initial state", box2.label = "1: progression",
  box3.label = "2: death", arrowLabelSymbol = "h"
)

## -----------------------------------------------------------------------------
library(simIDM)
transitionTrt <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
transitionPl <- exponential_transition(h01 = 0.5, h02 = 0.3, h12 = 0.6)

transitionList <- list(transitionPl, transitionTrt)

## -----------------------------------------------------------------------------
timepoints <- c(0, 0.1, 0.3, 0.7, 1, 5)
# OS Survival function for Constant transition hazards:
ExpSurvOS(timepoints, h01 = 0.2, h02 = 0.4, h12 = 0.1)
# OS Survival function for Weibull transition hazards:
WeibSurvOS(timepoints, h01 = 0.2, h02 = 0.5, h12 = 2.1, p01 = 1.2, p02 = 0.9, p12 = 1)
# OS Survival function for Piecewise Constant transition hazards:
PWCsurvOS(timepoints,
  h01 = c(0.3, 0.5), h02 = c(0.5, 0.8), h12 = c(0.7, 1),
  pw01 = c(0, 4), pw02 = c(0, 8), pw12 = c(0, 3)
)

## -----------------------------------------------------------------------------
timepoints <- c(0, 0.1, 0.3, 0.7, 1, 5)
# PFS Survival function for Constant transition hazards:
ExpSurvPFS(timepoints, h01 = 0.2, h02 = 0.4)
# PFS Survival function for Weibull transition hazards:
WeibSurvPFS(timepoints, h01 = 0.2, h02 = 0.5, p01 = 1.2, p02 = 0.9)
# PFS Survival function for Piecewise Constant transition hazards:
PWCsurvPFS(timepoints,
  h01 = c(0.3, 0.5), h02 = c(0.5, 0.8),
  pw01 = c(0, 4), pw02 = c(0, 8)
)

## -----------------------------------------------------------------------------
hTrtPFS <- sum(unlist(transitionTrt$hazards[1:2]))
hPlPFS <- sum(unlist(transitionPl$hazards[1:2]))
hRatioPFS <- hTrtPFS / hPlPFS
hRatioPFS

## -----------------------------------------------------------------------------
hRatioOS <- avgHRExpOS(transitionByArm = transitionList, alpha = 0.5, upper = 1000)
hRatioOS

## ----results = 'hide'---------------------------------------------------------
transitionListNull <- list(transitionPl, transitionPl)
nRep <- 100
simNull <- getClinicalTrials(
  nRep = nRep, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionListNull,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 100)
)

## -----------------------------------------------------------------------------
alphaOS <- 0.04
alphaPFS <- 0.01
criticalOS <- abs(qnorm(alphaOS / 2))
criticalPFS <- abs(qnorm(alphaPFS / 2))

## -----------------------------------------------------------------------------
library(rpact)
# Number of PFS events required for 80 % power via Schoenfeld:
schoenfeldPFS <- getSampleSizeSurvival(
  lambda2 = hPlPFS, hazardRatio = hRatioPFS,
  accrualTime = 0, accrualIntensity = 500,
  maxNumberOfSubjects = 1000, sided = 1,
  alpha = alphaPFS / 2, beta = 0.2
)
eventNumPFS <- ceiling(schoenfeldPFS$eventsFixed)
eventNumPFS

# Number of OS events required for 80 % power via Schoenfeld:
schoenfeldOS <- getSampleSizeSurvival(
  lambda2 = hPlPFS, hazardRatio = hRatioOS,
  accrualTime = 0, accrualIntensity = 500,
  maxNumberOfSubjects = 1000, sided = 1,
  alpha = alphaOS / 2, beta = 0.2
)
eventNumOS <- ceiling(schoenfeldOS$eventsFixed)
eventNumOS

## -----------------------------------------------------------------------------
empSignificant(
  simTrials = simNull,
  criticalPFS = criticalPFS,
  criticalOS = criticalOS,
  eventNumPFS = eventNumPFS,
  eventNumOS = eventNumOS
)

## ----results = 'hide'---------------------------------------------------------
simH1 <- getClinicalTrials(
  nRep = nRep, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionList,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 100)
)

## -----------------------------------------------------------------------------
empSignificant(
  simTrials = simH1,
  criticalPFS = criticalPFS,
  criticalOS = criticalOS,
  eventNumPFS = eventNumPFS,
  eventNumOS = eventNumOS
)

## -----------------------------------------------------------------------------
# median time point for 329 PFS events to have occurred:
timePointsPFS <- lapply(simH1, getTimePoint,
  eventNum = 329, typeEvent = "PFS",
  byArm = FALSE
)
median(unlist(timePointsPFS))

# median time point for 684 OS events to have occurred:
timePointsOS <- lapply(simH1, getTimePoint,
  eventNum = 684, typeEvent = "OS",
  byArm = FALSE
)
median(unlist(timePointsOS))

# mean number of PFS events at time of OS analysis
eventsPFS <- lapply(
  seq_along(timePointsPFS),
  function(t) {
    sum(simH1[[t]]$OSevent[(simH1[[t]]$OStime + simH1[[t]]$recruitTime) <= timePointsPFS[[t]]])
  }
)
mean(unlist(eventsPFS))

# mean number of OS events at time of PFS analysis
eventsOS <- lapply(
  seq_along(timePointsOS),
  function(t) {
    sum(simH1[[t]]$PFSevent[(simH1[[t]]$PFStime + simH1[[t]]$recruitTime) <= timePointsOS[[t]]])
  }
)
mean(unlist(eventsOS))

