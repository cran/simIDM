## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE,  fig.cap = "Figure 1 - Multistate model with indermediate state progession and absorbing state death",fig.height = 5, fig.width = 8, out.width = "60%"----
library(prodlim)
plotIllnessDeathModel(
  style = 1, box1.label = "0: initial state", box2.label = "1: progression",
  box3.label = "2: death", arrowLabelSymbol = "h"
)

## -----------------------------------------------------------------------------
nPat <- c(30, 60)

## -----------------------------------------------------------------------------
library(simIDM)

transitionGroup1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
print(transitionGroup1$hazards)
transitionGroup2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
print(transitionGroup2$hazards)
transitionbyArm <- list(transitionGroup1, transitionGroup2)

## -----------------------------------------------------------------------------
dropout <- list(rate = 0.1, time = 12)

## -----------------------------------------------------------------------------
# first example
accrual <- list(param = "intensity", value = 12)
# second example
accrual <- list(param = "intensity", value = 3)

## -----------------------------------------------------------------------------
simStudies1 <- getClinicalTrials(
  nRep = 100, nPat = nPat, seed = 1234, datType = "1rowTransition",
  transitionByArm = transitionbyArm, dropout = dropout,
  accrual = accrual
)

simStudies2 <- getClinicalTrials(
  nRep = 100, nPat = nPat, seed = 1234, datType = "1rowPatient",
  transitionByArm = transitionbyArm, dropout = dropout,
  accrual = accrual
)

## -----------------------------------------------------------------------------
head(simStudies1[[1]], 6)
head(simStudies2[[1]], 5)

## -----------------------------------------------------------------------------
simStudy_Trans <- getOneClinicalTrial(
  nPat = nPat, transitionByArm = transitionbyArm,
  dropout = dropout,
  accrual = accrual
)

simStudy_PFSOS <- getDatasetWideFormat(simStudy_Trans)
head(simStudy_Trans, 10)
head(simStudy_PFSOS, 10)

## -----------------------------------------------------------------------------
simStudy_Trans <- getOneClinicalTrial(
  nPat = nPat, transitionByArm = transitionbyArm,
  dropout = dropout,
  accrual = accrual
)

simStudy_PFSOS <- getDatasetWideFormat(simStudy_Trans)
head(simStudy_Trans, 10)
head(simStudy_PFSOS, 10)

simStudy_40PFS <- censoringByNumberEvents(data = simStudy_PFSOS, eventNum = 40, typeEvent = "PFS")
simStudy_30POS <- censoringByNumberEvents(data = simStudy_PFSOS, eventNum = 30, typeEvent = "PFS")
simStudy_40PFS
simStudy_30POS

