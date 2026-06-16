# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

Bartlett_fn <- function(x) ifelse(abs(x)>=1,0,1-abs(x))

