# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

cov_est_fn <- function(tau,Y,X,Z,Lambda,Lambda.derivative,beta.hat,Itilde,Itilde.deriv,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE,h.adj=1) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  G.hat <- G_est_fn(Y=Y,X=X,Z=Z,Lambda=Lambda,Lambda.derivative=Lambda.derivative,beta.hat=beta.hat,Itilde.deriv=Itilde.deriv,h=h^h.adj,VERBOSE=VERBOSE)
  Ginv <- tryCatch(solve(G.hat), error=function(w)NA)
  if (is.na(Ginv[1])) return(NA) else {
    LRV.hat <- LRV_est_fn(tau=tau,Y=Y,X=X,Z=Z,Lambda=Lambda,beta.hat=beta.hat,Itilde=Itilde,h=h,structure=structure,cluster.X.col=cluster.X.col,LRV.kernel=LRV.kernel,LRV.ST=LRV.ST,VERBOSE=VERBOSE)
    return(Ginv %*% LRV.hat %*% t(Ginv))
  }
}

