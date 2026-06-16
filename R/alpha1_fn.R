# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

alpha1_fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/((1-rhos)^6*(1+rhos)^2)) / sum(ws*sigmas^4/(1-rhos)^4)

