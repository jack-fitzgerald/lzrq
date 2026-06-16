# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

Itilde_KS17 <- function(u) {  ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0)) }

