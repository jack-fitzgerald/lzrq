# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

gmmq_wald_test <- gmmq_Wald_test <- function(a.hat,A.hat,tau,Y,X,Z,Lambda,Lambda.derivative,beta.hat,Itilde,Itilde.deriv,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE,h.adj=1) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  if (any(is.na(Z))) stop("argument Z should not have any NA entries. Use ret <- gmmq(...,RETURN.Z=TRUE); gmmq_wald_test(...,Z=ret$Z,...)")
  NAany <- which(apply(X=cbind(Y,X),MARGIN=1,FUN=function(row)any(is.na(row))))
  if (length(NAany)>0) {
    cat(sprintf("Removing %d observations due to NA in Y or X\n",length(NAany)))
    Y <- matrix(Y[-NAany,], ncol=dim(Y)[2])
    X <- matrix(X[-NAany,], ncol=dim(X)[2])
  }
  # 
  n <- dim(Z)[1]
  if (dim(Y)[1]!=n || dim(X)[1]!=n) stop("The number of non-NA rows in Y, X, and Z should be the same.")
  # 
  cov.hat <- cov_est_fn(tau=tau,Y=Y,X=X,Z=Z,Lambda=Lambda,Lambda.derivative=Lambda.derivative,beta.hat=beta.hat,Itilde=Itilde,Itilde.deriv=Itilde.deriv,h=h,structure=structure,cluster.X.col=cluster.X.col,LRV.kernel=LRV.kernel,LRV.ST=LRV.ST,VERBOSE=VERBOSE,h.adj=h.adj)
  if (is.na(cov.hat[1])) {
    return(data.frame(Wald.stat=NA, pval=NA, df=length(a.hat)))
  } else {
    W.hat <- tryCatch(n * matrix(a.hat,nrow=1) %*% solve(matrix(A.hat,nrow=length(a.hat)) %*% cov.hat %*% t(matrix(A.hat,nrow=length(a.hat)))) %*% matrix(a.hat,ncol=1), error=function(w)NA)
    if (is.na(W.hat[1])) return(data.frame(Wald.stat=NA, pval=NA, df=length(a.hat)))
    r <- length(a.hat)
    pval <- 1 - pchisq(q=W.hat, df=r)
    return(data.frame(Wald.stat=W.hat,pval=pval,df=r))
  }
}

