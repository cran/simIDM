## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(simIDM)

# constant hazards:
transitionExp <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)

# Weibull hazards:
transitionWeib <- weibull_transition(h01 = 1, h02 = 1.2, h12 = 1.3, p01 = 1.1, p02 = 0.8, p12 = 1.2)

# piecewise constant hazards:
transitionPwc <- piecewise_exponential(
  h01 = c(1, 1.3), h02 = c(0.8, 1.5), h12 = c(1, 1),
  pw01 = c(0, 3), pw02 = c(0, 1), pw12 = c(0, 8)
)

## -----------------------------------------------------------------------------
# constant hazards:
corTrans(transitionExp)

# Weibull hazards:
corTrans(transitionWeib)

# piecewise constant hazards:
corTrans(transitionPwc)

## -----------------------------------------------------------------------------
transitionExp <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
simData <- getOneClinicalTrial(
  nPat = c(500), transitionByArm = list(transitionExp),
  dropout = list(rate = 0.8, time = 12),
  accrual = list(param = "time", value = 1)
)

## -----------------------------------------------------------------------------
# Create TransitionParameters object with starting values for ML estimation:
transition <- exponential_transition(h01 = 1, h02 = 1, h12 = 1)
# Estimate parameters:
est <- estimateParams(data = simData, transition = transition)
# Get estimated transition hazards:
est$hazards

## -----------------------------------------------------------------------------
corPFSOS(data = simData, transition = transition, bootstrap = TRUE, conf_level = 0.95)

