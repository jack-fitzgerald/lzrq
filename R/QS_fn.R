# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

QS_fn <- function(x) ifelse(x==0,1,(25/(12*pi^2*x^2))*(sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5)))

