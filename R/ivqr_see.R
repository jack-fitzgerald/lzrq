# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

ivqr_see <- function(tau, Y, D=NULL, X.exog=NULL, Z.excl=NULL, h, b.init=NULL, iidSE=TRUE) {
  # Validate/set up data
  dD <- 0
  if (!is.null(D)) {
    D <- as.matrix(D)
    dD <- ncol(D)
  }
  if (is.null(X.exog)) {
    X <- matrix(data=1, nrow=length(Y), ncol=1)
  } else {
    X.exog <- as.matrix(X.exog)
    X <- cbind(X.exog,1)
  }
  dX <- ncol(X)
  if (is.null(Y)) stop('Must specify Y argument.') else Y <- as.matrix(Y)
  if (!is.null(Z.excl)) Z.excl <- as.matrix(Z.excl)
  # Set b.init (if not specified by user): QR if available, else 0
  if (is.null(b.init)) {
    if (require(quantreg)) {
      if (is.null(X.exog)) {
        tmp <- coef(rq(Y~D))
      } else if (is.null(D)) {
        tmp <- coef(rq(Y~X.exog))
      } else {
        tmp <- coef(rq(Y~D+X.exog))
      }
      b.init <- tmp[c(2:length(tmp),1)]
      if (!is.numeric(b.init) || any(is.nan(b.init)) || any(is.infinite(b.init)) || length(b.init)<1) b.init <- 0
    } else {
      b.init <- 0
    }
  }
  # Set bandwidth (if not specified by user)
  if (missing(h)) {
    ret0 <- gmmq(tau=tau, Y=cbind(Y,D), X=X, 
                 Z.excl=Z.excl, dB=dD+dX, b.init=b.init)
    if (length(b.init)==1 && b.init==0) b.init <- ret0$b
    hhat <- ivqr_bw(p=tau, Y=Y, X=cbind(D,X), 
                    b.init=ret0$b)
    if (!is.numeric(hhat) || is.nan(hhat) || is.infinite(hhat)) {
      warning("Problem with plug-in bandwidth; instead using smallest feasible bandwidth.")
      return(list(b=ret0$b, h=ret0$h, hhat=ret0$h))
    }
  } else {
    hhat <- h
  }
  # Actually run the estimator
  ret1 <- gmmq(tau=tau, Y=cbind(Y,D), X=X, 
               Z.excl=Z.excl, dB=dD+dX, h=hhat, b.init=b.init, iidSE=iidSE)
  # Add names/labels for coefficients
  ret1$b <- c(ret1$b)
  if (is.null(D)) {
    names(ret1$b) <- c(sprintf('exog.%d',1:ncol(cbind(X.exog))),'(Intercept)')
  } else if (is.null(X.exog)) {
    names(ret1$b) <- c(sprintf('endog.%d',1:ncol(cbind(D))),'(Intercept)')
  } else {
    names(ret1$b) <- c(sprintf('endog.%d',1:ncol(cbind(D))),
                       sprintf('exog.%d',1:ncol(cbind(X.exog))),'(Intercept)')
  }
  ret <- list(b=ret1$b, h=ret1$h, hhat=hhat)
  if (iidSE) ret$se <- ret1$se
  return(ret)
}

