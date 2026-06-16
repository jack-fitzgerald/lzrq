# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

LRV_est_fn <- function(tau,Y,X,Z,Lambda,beta.hat,Itilde,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  n <- dim(Z)[1]
  if (structure %in% c('iid','ts')) {
    if (structure=='iid') {
      LRV.kernel <- 'uniform'; LRV.lag <- 0; LRV.ST <- 1
    }
    if (missing(LRV.kernel) || !is.character(LRV.kernel)) stop("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS' when structure is 'ts'")
    if (LRV.kernel=='uniform') weight.fn <- uniform_fn else if (LRV.kernel=='Bartlett') weight.fn <- Bartlett_fn else if (LRV.kernel=='QS') weight.fn <- QS_fn else stop(sprintf("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS'; not %s",LRV.kernel))
    # Compute gni() matrix
    gni.mat <- Z*array(data=Itilde(-Lambda(y=Y,x=X,b=beta.hat)/h)-tau,dim=dim(Z))
    #
    if (is.na(LRV.ST)) { # Set ST automatically
      rho.hats <- sigma.hats <- rep(NA,dim(Z)[2])
      for (a in 1:length(rho.hats)) {
        rho.hats[a] <- sum(gni.mat[1:(n-1),a]*gni.mat[2:n,a]) / sum(gni.mat[1:(n-1),a]^2)
        sigma.hats[a] <- suppressWarnings(sqrt(var(gni.mat[,a]) * (1-rho.hats[a]^2)))
      }
      if (any(c(is.nan(c(rho.hats,sigma.hats)),is.na(c(rho.hats,sigma.hats))))) {
        LRV.ST <- n^(1/5); if (LRV.kernel=='uniform' || LRV.kernel=='Bartlett') LRV.ST <- n^(1/3)
        warning(sprintf("AR method from Andrews (1991) for selecting S_T returned NA or NaN values; using S_T=%g",LRV.ST))
      } else if (LRV.kernel=='uniform') {
        LRV.ST <- uniform_ST_fn(alpha2=alpha2_fn(rhos=rho.hats,sigmas=sigma.hats),n=n)
      } else if (LRV.kernel=='Bartlett') {
        LRV.ST <- Bartlett_ST_fn(alpha1=alpha1_fn(rhos=rho.hats,sigmas=sigma.hats),n)
      } else if (LRV.kernel=='QS') {
        LRV.ST <- QS_ST_fn(alpha2=alpha2_fn(rhos=rho.hats,sigmas=sigma.hats),n)
      } else stop("Uncaught case.")
    }
    if (LRV.kernel!='QS') LRV.lag <- floor(LRV.ST)
    #
    tmpsum <- array(0,dim=rep(dim(Z)[2],2))
    for (i in 1:n) {
      if (LRV.kernel=='QS') krange <- 1:n else krange <- max(1,i-LRV.lag):min(n,i+LRV.lag)
      for (k in krange) {
        tmpsum <- tmpsum + 
          weight.fn((i-k)/LRV.ST) * 
          (matrix(gni.mat[i,],ncol=1) %*% matrix(gni.mat[k,],nrow=1))
      }
    }
    return(tmpsum/(n-length(beta.hat))) #denominator adjustment per Andrews (1991) eqn (2.5)
  } else if (structure=='cluster') {
    stop("Not yet implemented: clustered covariance estimation")
  } else stop(sprintf("Argument structure must be either 'iid' or 'ts' or 'cluster' but its value is %s",structure))
}

