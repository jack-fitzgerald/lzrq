# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

gmmq_BS_test <- function(a.fn,A.fn,tau,Y,X,Z,Z.excl,Lambda,Lambda.derivative,beta.hat,beta0, Itilde=NULL, Itilde.deriv=NULL, h,structure=c('iid','ts','cluster'),cluster.X.col=NULL,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE,h.adj=1,BREP=99,BLAG=NA) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  if (is.null(Itilde)) Itilde <- get("Itilde_KS17", envir = parent.frame())
  if (is.null(Itilde.deriv)) Itilde.deriv <- get("Itilde_deriv_KS17", envir = parent.frame())
  if (any(is.na(Z))) stop("argument Z should not have any NA entries. Use ret <- gmmq(...,RETURN.Z=TRUE); gmmq_wald_test(...,Z=ret$Z,...)")
  NAany <- which(apply(X=cbind(Y,X,Z.excl),MARGIN=1,FUN=function(row)any(is.na(row))))
  if (length(NAany)>0) {
    cat(sprintf("Removing %d observations due to NA in Y or X\n",length(NAany)))
    Y <- matrix(Y[-NAany,], ncol=dim(Y)[2])
    X <- matrix(X[-NAany,], ncol=dim(X)[2])
    Z.excl <- matrix(Z.excl[-NAany,], ncol=dim(Z.excl)[2])
  }
  # 
  # Set RNG seed for replicability (and save current seed to re-seed on exit)
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  set.seed(112358) #for replicability
  # 
  n <- dim(Z)[1]
  if (dim(Y)[1]!=n || dim(X)[1]!=n || dim(Z.excl)[1]!=n) stop("The number of non-NA rows in Y, X, Z.excl, and Z should be the same.")
  # 
  W.hat <- gmmq_wald_test(a.hat=a.fn(beta.hat,beta0),A.hat=A.fn(beta.hat,beta0),tau=tau,Y=Y,X=X,Z=Z,Lambda=Lambda,Lambda.derivative=Lambda.derivative,beta.hat=beta.hat,Itilde=Itilde,Itilde.deriv=Itilde.deriv,h=h,structure=structure,cluster.X.col=cluster.X.col,LRV.kernel=LRV.kernel,LRV.ST=LRV.ST,VERBOSE=VERBOSE,h.adj=h.adj)$Wald.stat
  W.hat.stars <- rep(NA,BREP)
  if (structure=='cluster') {
    cluster.vals <- unique(X[,cluster.X.col])
    n.cluster <- length(cluster.vals)
    clustinds <- vector("list",n.cluster)
    for (i in 1:n.cluster) clustinds[[i]] <- which(X[,cluster.X.col]==cluster.vals[i])
  }
  for (b in 1:BREP) {
    if (structure=='iid') {
      indstars <- sample(1:n,n,TRUE)
    } else if (structure=='ts') {
      # stationary block bootstrap (Politis and Romano, 1994): circular, but random block lengths
      BS.p <- n^(-1/3) #see p. 1306 of Politis and Romano (1994)
      if (!is.na(BLAG)) BS.p <- 1/BLAG
      indstars <- rep(NA,n)
      indstars[1] <- sample(1:n,1)
      for (i in 2:n) {
        indstars[i] <- ifelse(runif(1)<BS.p,sample(1:n,1),ifelse(indstars[i-1]<n,indstars[i-1]+1,1)) # p. 1304
      }
    } else if (structure=='cluster') {
      indstars <- unlist(clustinds[sample(1:n.cluster,n.cluster,TRUE)])
    } else stop(sprintf("Argument 'structure' to function gmmq_BS_test() must be 'iid' or 'ts' or 'cluster'"))
    Ystar <- matrix(Y[indstars,],ncol=ncol(Y))
    Xstar <- matrix(X[indstars,],ncol=ncol(X))
    Zexclstar <- matrix(Z.excl[indstars,],ncol=ncol(Z.excl))
    # 
    retstar <- suppressWarnings(tryCatch(gmmq(tau=tau,Y=Ystar,X=Xstar,Z.excl=Zexclstar,dB=length(beta.hat),
                                              Lambda=Lambda, Lambda.derivative=Lambda.derivative,
                                              h=h, VERBOSE=VERBOSE, RETURN.Z=TRUE,
                                              b.init=beta.hat),
                                         error=function(w)NA) )
    if (!is.na(retstar)[1]) W.hat.stars[b] <- 
      suppressWarnings(gmmq_wald_test(a.hat=a.fn(retstar$b,beta.hat), A.hat=A.fn(retstar$b,beta.hat), tau=tau, Y=Ystar, X=Xstar, Z=retstar$Z, Lambda=Lambda,Lambda.derivative=Lambda.derivative, beta.hat=retstar$b, Itilde=Itilde,Itilde.deriv=Itilde.deriv, h=retstar$h, structure=structure, cluster.X.col=cluster.X.col, LRV.kernel=LRV.kernel, LRV.ST=LRV.ST, VERBOSE=VERBOSE, h.adj=h.adj)$Wald.stat)
  }
  pval.BS <- mean(W.hat.stars>W.hat,na.rm=TRUE)
  if (sum(is.na(W.hat.stars))>0) warning(sprintf("Some NA replications; only %d reps being used. (Increase BREP to increase #non-NA reps.)",sum(!is.na(W.hat.stars))))
  return(list(pval=pval.BS, Wald.stat=W.hat, BS.Wald.stats=W.hat.stars))
}

