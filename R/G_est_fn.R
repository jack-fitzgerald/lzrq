# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

G_est_fn <- function(Y,X,Z,Lambda,Lambda.derivative,beta.hat,Itilde.deriv,h,VERBOSE=FALSE) {
  n <- dim(Z)[1]
  L <- Lambda(Y,X,beta.hat)
  Ld <- Lambda.derivative(Y,X,beta.hat)
  tmpsum2 <- t(array(data=Itilde.deriv(-L/h),dim=dim(Z)) * Z) %*% Ld
  return(-tmpsum2/(n*h))
}

