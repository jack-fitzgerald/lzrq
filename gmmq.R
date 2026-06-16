# "Smoothed GMM for quantile models" by de Castro, Galvao, Kaplan, and Liu
# Journal of Econometrics, 2019: https://doi.org/10.1016/j.jeconom.2019.04.008
# This code covers exact identification; see also ivqr_gmm.R
# Questions? Comments?  contact Dave (https://kaplandm.github.io)

# pracma::newtonsys used for solving sample moment conditions (Z-estimator)
if (!require(pracma)) stop('Please install (from CRAN) the R package pracma to run this code: install.packages("pracma")')

# General idea: estimator based on moment conditions E{z*[1{Lambda(y,x,b)<=0}-tau]}=0
# tau: quantile index between zero and one, e.g. tau=0.5 for median.
# dB: length of vector of unknown parameters.
# Y: n-by-dY matrix of endogenous variables, where n=sample size.
# X: n-by-dX matrix of exogenous variables used in function Lambda. Note: just use X=matrix(1,ncol=1,nrow=n) if not really needed. Note: X need not include a constant since the function Lambda can add a constant, *but* in the case of overidentification it may be helpful to explicitly include a column of 1's as the first column in X, as well as any transformations of other columns that are used in Lambda; see note on overidentification below.
# Z.excl: if applicable, n-by-dZ matrix of excluded instruments (i.e., exogenous but not in Lambda); else, use NULL (the default). If dZ=dB, then only Z.excl (not X) is used to construct the moment conditions.
# NOTES on identification: 
# must have dB<=dX+dZ (necessary [but not sufficient] condition for the parameter vector to be identified)
# Overidentification: if dB<dX+dZ (strictly), then moment conditions are transformed until the number of parameters (dB) equals the number of moments. If dB equals the number of columns in Y plus the number of columns in X minus one, i.e. dB==dim(Y)[2]+dim(X)[2]-1, then we take a linear projection of cbind(Y[,-1],X) onto cbind(X,Z.excl) to get the full set of dB instruments. For example, if Y has two columns, Y=cbind(Y1,Y2) where Y1 is an outcome and Y2 is an endogenous regressor, and X contains exogenous regressors (including a constant, and any transformations used in Lambda), and Z.excl is n-by-dZ with dZ>1 (so there is overidentification), then we project cbind(Y2,X) onto cbind(X,Z.excl); see Kaplan and Sun (2017, Econometric Theory), "Smoothed estimating equations for instrumental variables quantile regression."  If this condition dB==dim(Y)[2]+dim(X)[2]-1 does not hold, then some of the Z.excl are simply removed until identification is exact; this is not optimal, but still yields a consistent estimator. Ideally, the user will determine the best set of instruments in an exactly identified system and adjust Z.excl accordingly; if dZ=dB, then only Z.excl (not X) is used for the unconditional moments.
# apparent lack of identification: if you have dB>dX+dZ, but you have a *conditional* moment condition (or, conditional quantile restriction), then you should add transformations of the original exogenous variables/instruments until dB=dX+dZ. Example: y1 is the outcome, y2 is the endogenous regressor, and z.excl is an instrument; all are scalars. A random coefficients structural model y1=a+(y2,y2^2)*(b1(U),b2(U)) with unobserved U leads (under some assumptions) to the conditional quantile restriction P{y1<a+y2*b1(tau)+y2^2*b2(tau) | z}=tau as in Theorem 1 of Chernozhukov and Hansen (2005). In terms of our code, we'd initially have X=matrix(1,nrow=n,ncol=1), an n-by-1 array of 1's; Y=cbind(Y1,Y2,Y2^2), an n-by-3 matrix; and Z.excl=Z, an n-by-1 array also. (Or, squaring y2 could be done inside Lambda.) However, then dB=3, dX=1, and dZ=1, so dB>dX+dZ since 3>2. However, since we have a *conditional* (on z) quantile condition, we can add a transformation of z to get another *unconditional* moment condition and then dB=dX+dZ, e.g. setting Z.excl=cbind(Z,Z^2) so that dZ=2 and thus 3=1+2.
# Lambda: a function of matrices y and x and vector b, where dim(y)=(n,dY), dim(x)=(n,dX), and b has length dB. The moment condition is E{z*[1{Lambda(y,x,b)<=0}-tau]}=0 where vector z includes both x and z.excl.  Example: Lambda <- function(y,x,b) y-(x%*%b)  or y[,1]-cbind(y[,-1],x)%*%b
# Lambda.derivative: a function of matrices y and x and vector b that returns an n-by-dB matrix of (partial) first derivatives of Lambda with respect to b, evaluated at each observation.  This greatly help the numerical solver and is needed for hypothesis testing.  Example: if Lambda(y,x,b)=y-x%*%b, then Lambda.derivative=function(y,x,b) { -x }   Can be NULL.
# h: initial smoothing bandwidth to try.  If h=0 (default), then the code tries to find the smallest possible bandwidth.
# VERBOSE: output update messages each iteration of the search over h
# RETURN.Z: if TRUE then also return the used instrument matrix Z; may be useful for subsequently calling gmmq.wald.test(), but otherwise set to FALSE
# b.init: starting value for the vector beta in the numerical search.
# 
# Return value: list with estimate of the b parameter vector and actual h used (and the instrument matrix Z, if RETURN.Z is TRUE)
# 
# Smoothing function Itilde() and its derviative Itilde'(); Vector input supported
# Suggestion from Kaplan and Sun (2017):
Itilde.KS17 <- function(u) {  ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0)) }
Itilde.deriv.KS17 <- function(u) { ifelse(u > -1 & u < 1, (105/64)*(1-5*u^2+7*u^4-3*u^6), 0) }
#
gmmq <- function(tau, Y, X, Z.excl=NULL, dB, 
                 Lambda=function(y,x,b) y[,1]-cbind(y[,-1],x)%*%b, 
                 Lambda.derivative=function(y,x,b) -cbind(y[,-1],x), 
                 h=0, VERBOSE=FALSE, RETURN.Z=FALSE, b.init=0, iidSE=FALSE) {
  # pracma::newtonsys used for solving sample moment conditions (Z-estimator)
  if (!require(pracma)) stop('Please install (from CRAN) the R package pracma to run this code.')
  
  # Check arguments: errors, warnings, defaults
  if (missing(tau) || missing(Y) || missing(X)) stop('Must supply tau (scalar), dB (scalar), Y (matrix), X (matrix)')
  if (!is.numeric(tau) || !is.numeric(h)) stop('Arguments tau and h must be numeric.')
  Y <- as.matrix(Y);  X <- as.matrix(X)
  if (!is.null(Z.excl)) Z.excl <- as.matrix(Z.excl)
  if (missing(dB)) dB <- dim(Y)[2]+dim(X)[2]-1 # for linear IVQR
  if (!is.numeric(dB)) stop("dB must be numeric (or left unspecified, in which case it defaults to dim(Y)[2]+dim(X)[2]-1 for linear IVQR)")
  if (b.init[1]!=0 && length(b.init)!=dB) stop(sprintf("Length of b.init must equal dB, which is %d", dB))
  if (!is.function(Lambda)) stop('Lambda must be a function.')
  if (!is.null(Lambda.derivative) && !is.function(Lambda.derivative)) stop('Lambda.derivative must be a function (or NULL).')
  if (iidSE && is.null(Lambda.derivative)) stop('Must set iidSE=FALSE if Lambda.derivative=NULL')
  if (!is.array(Y) || !is.array(X) || (!is.null(Z.excl) && !is.array(Z.excl))) stop('Y, X, and Z.excl (if not NULL) must be arrays/matrices, so is.array(Y)=TRUE, etc. To make a vector into a column array, use e.g. matrix(Y,ncol=1)')
  if (length(tau)!=1 || length(dB)!=1 || length(h)!=1) stop('tau, dB, and h must all be scalars: length(dB)==1, etc.')
  dY <- dim(Y)[2]
  if (is.null(Z.excl)) Z <- X else {
    if (dim(X)[1]!=dim(Z.excl)[1]) stop('X and Z.excl must have same number of rows.') else Z <- cbind(X,Z.excl)
  }
  if (dim(Z)[1]!=dim(Y)[1]) stop('Y and X and Z.excl must have same number of rows.')
  if (tau<=0 || tau>=1) stop('Quantile index tau (first argument) must be strictly between 0 and 1.')
  if (h<0) stop('Initial bandwidth h must be nonnegative.')
  # Check for NA, remove rows
  NAany <- which(apply(X=cbind(Y,Z),MARGIN=1,FUN=function(row)any(is.na(row))))
  if (length(NAany)>0) {
    cat(sprintf("Removing %d observations due to NA in Y, X, or Z.excl\n",length(NAany)))
    Y <- matrix(Y[-NAany,], ncol=dim(Y)[2])
    Z <- matrix(Z[-NAany,], ncol=dim(Z)[2])
    X <- matrix(X[-NAany,], ncol=dim(X)[2])
    Z.excl <- matrix(Z.excl[-NAany,], ncol=dim(Z.excl)[2])
  }
  n <- dim(Y)[1]
  #
  if (dB>dim(Z)[2]) stop('Parameters not identified: must have dB<=dim(X)[2]+dim(Z.excl)[2]') else if (dB<dim(Z)[2]) { #overidentified
    if (dB==dim(Y)[2]+dim(X)[2]-1) { # for consistency with IVQR-SEE
      tmp1 <- lm(Y[,-1]~Z) #could Z+0 but just in case Z excludes intercept...
      tmp2 <- lm(X~Z)
      Z <- cbind(predict(tmp1),predict(tmp2))
    } else if (dim(Z.excl)[2]==dB) { # user supplied all dB instruments in Z.excl
      Z <- Z.excl
    } else { #for now, just ignore some Z.excl
      Z <- matrix(Z[,1:dB], ncol=dB)
    }
  }
  
  # Set smoothing function Itilde() and its derviative Itilde'()
  # (Vector input supported)
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
  
  # For solver (to find solution to estimating equations)
  obj.fn <- function(b,h) { #objective function
    L <- Lambda(y=Y,x=X,b=b)
    return(colMeans(Z*array(data=Itilde(-L/h)-tau,dim=dim(Z))))
  }
  if (!is.null(Lambda.derivative)) {
    jac.fn <- function(b,h) { #Jacobian function
      L <- Lambda(y=Y,x=X,b=b)
      Lp <- Lambda.derivative(y=Y,x=X,b=b)
      return((t(Z) %*% (Lp*array(data=Itilde.deriv(-L/h),dim=dim(Z))*(-1/h)))/n)
    }
  }
  
  # Compute estimator: if h=0, find smallest h newtonsys can use; else use user's h (if big enough)
  if (h==0) {
    h.init <- h.cur <- 0.001 #was: .Machine$double.eps * 1e3
    hfac <- 100
  } else {
    h.init <- h.cur <- h
    hfac <- 10
  }
  h.lo <- h
  h.hi <- NA
  b.best <- matrix(data=b.init,nrow=dB,ncol=1)
  MAXITER <- 400;  H.MAX <- 1e10 #.Machine$double.xmax/1e3
  # Loop: h<=h.lo<h.cur<h.hi, where h.hi can be solved by newtonsys but h.lo cannot
  if (VERBOSE) cat(sprintf('n=%d, tau=%g, dB=%d, desired h=%g\n',n,tau,dB,h))
  last.error <- NULL
  while (h.cur<=H.MAX && h.cur>.Machine$double.eps*1e2 && 
         (is.na(h.hi) || (h.hi/h.lo)>1.4)) {
    if (VERBOSE) cat(sprintf("h.cur=%g: ",h.cur))
    if (is.null(Lambda.derivative)) J <- NULL else J <- function(b)jac.fn(b,h.cur)
    soln <- tryCatch(newtonsys(Ffun=function(b)obj.fn(b,h.cur),Jfun=J,
                               x0=b.best,maxiter=MAXITER),
                     warning=function(w) NA,
                     error=function(e) { list(NA,e) } #warning("(above) error from running pracma::newtonsys."); stop(e) } #NA
    )
    exitOK <- !(is.na(soln[1]) || soln$niter==MAXITER)
    if (!exitOK) {
      last.error <- soln[2]
      if (VERBOSE) {
        if (is.na(soln[1])) cat(sprintf("no solution, soln=NA\n")) else cat(sprintf("no solution, soln$niter=%d\n",soln$niter))
      }
      h.lo <- h.cur
      if (is.na(h.hi)) { #keep going up till find h.cur that actually works
        h.cur <- h.cur*hfac
      } else { #try arithmetic [not geometric] mean of h.lo and h.hi
        h.cur <- (h.lo+h.hi)/2 #sqrt(h.lo*h.hi)
      }
    } else { #newtonsys found a solution
      if (VERBOSE) cat(sprintf("Solution found!\n"))
      h.hi <- h.cur; b.best <- soln$zero; h.lo <- h
      if (h.cur==h) { #found solution w/ user-requested h
        break
      } else if (h.hi/h.init>2.5 || h==0) { 
        h.cur <- (h.lo+h.hi)*(2/3)
      } else h.cur <- h.init
    } 
  }
  if (h.cur>H.MAX) {
    warning("(above) error from pracma::newtonsys")
    stop(last.error)
  }
  if (iidSE) {
    L <- Lambda(y=Y, x=X, b=b.best)
    h.SE <- bw.nrd(L)
    J.hat <- jac.fn(b=b.best, h=h.SE)
    S.hat <- tau*(1-tau)* t(Z) %*% Z / n
    S.hat.inv <- solve(S.hat)
    cov.hat <- solve(t(J.hat) %*% S.hat.inv %*% J.hat) / n
  }
  # 
  ret <- list(b=b.best, h=h.hi)
  if (RETURN.Z) ret$Z <- Z
  if (iidSE) ret$se <- sqrt(diag(cov.hat))
  return(ret)
}





# Wald test of H0:a(beta)=0
# a.hat: the vector a(beta.hat)
# A.hat: the partial derivative (i.e. Jacobian) matrix of a() evaluated at beta.hat, where row i, column k entry is the partial derivative of the i-th element of a() wrt the k-th element of beta.
gmmq.wald.test <- gmmq.Wald.test <- function(a.hat,A.hat,tau,Y,X,Z,Lambda,Lambda.derivative,beta.hat,Itilde,Itilde.deriv,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE,h.adj=1) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  if (any(is.na(Z))) stop("argument Z should not have any NA entries. Use ret <- gmmq(...,RETURN.Z=TRUE); gmmq.wald.test(...,Z=ret$Z,...)")
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
  cov.hat <- cov.est.fn(tau=tau,Y=Y,X=X,Z=Z,Lambda=Lambda,Lambda.derivative=Lambda.derivative,beta.hat=beta.hat,Itilde=Itilde,Itilde.deriv=Itilde.deriv,h=h,structure=structure,cluster.X.col=cluster.X.col,LRV.kernel=LRV.kernel,LRV.ST=LRV.ST,VERBOSE=VERBOSE,h.adj=h.adj)
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

# Studentized bootstrap (see gmmq.BS2.test for non-studentized); BREP is number of bootstrap replications (more is better but slower)
# Note: now a.fn and A.fn are functions (of beta.hat and beta0), not just values like a.hat and A.hat
# Note: both Z and Z.excl are arguments; use RETURN.Z=TRUE in gmmq() for the former
# Note: cluster BS implementation is simple; may not work well for {few/small/unbalanced/etc.} clusters
# BLAG: how many lags to include in bootstrap block (although this is translated to 1/p for a stationary bootstrap; see p. 1310 in Politis and Romano, 1994)
gmmq.BS.test <- function(a.fn,A.fn,tau,Y,X,Z,Z.excl,Lambda,Lambda.derivative,beta.hat,beta0, Itilde=NULL, Itilde.deriv=NULL, h,structure=c('iid','ts','cluster'),cluster.X.col=NULL,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE,h.adj=1,BREP=99,BLAG=NA) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  if (is.null(Itilde)) Itilde <- get("Itilde.KS17", envir = parent.frame())
  if (is.null(Itilde.deriv)) Itilde.deriv <- get("Itilde.deriv.KS17", envir = parent.frame())
  if (any(is.na(Z))) stop("argument Z should not have any NA entries. Use ret <- gmmq(...,RETURN.Z=TRUE); gmmq.wald.test(...,Z=ret$Z,...)")
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
  W.hat <- gmmq.wald.test(a.hat=a.fn(beta.hat,beta0),A.hat=A.fn(beta.hat,beta0),tau=tau,Y=Y,X=X,Z=Z,Lambda=Lambda,Lambda.derivative=Lambda.derivative,beta.hat=beta.hat,Itilde=Itilde,Itilde.deriv=Itilde.deriv,h=h,structure=structure,cluster.X.col=cluster.X.col,LRV.kernel=LRV.kernel,LRV.ST=LRV.ST,VERBOSE=VERBOSE,h.adj=h.adj)$Wald.stat
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
    } else stop(sprintf("Argument 'structure' to function gmmq.BS.test() must be 'iid' or 'ts' or 'cluster'"))
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
      suppressWarnings(gmmq.wald.test(a.hat=a.fn(retstar$b,beta.hat), A.hat=A.fn(retstar$b,beta.hat), tau=tau, Y=Ystar, X=Xstar, Z=retstar$Z, Lambda=Lambda,Lambda.derivative=Lambda.derivative, beta.hat=retstar$b, Itilde=Itilde,Itilde.deriv=Itilde.deriv, h=retstar$h, structure=structure, cluster.X.col=cluster.X.col, LRV.kernel=LRV.kernel, LRV.ST=LRV.ST, VERBOSE=VERBOSE, h.adj=h.adj)$Wald.stat)
  }
  pval.BS <- mean(W.hat.stars>W.hat,na.rm=TRUE)
  if (sum(is.na(W.hat.stars))>0) warning(sprintf("Some NA replications; only %d reps being used. (Increase BREP to increase #non-NA reps.)",sum(!is.na(W.hat.stars))))
  return(list(pval=pval.BS, Wald.stat=W.hat, BS.Wald.stats=W.hat.stars))
}



# non-Studentized bootstrap: return 2-sided symmetric CI for T.fn(beta), a scalar-valued function of the parameter vector; e.g., T.fn <- function(b) b[1] is used for a CI for the first element of the parameter vector.
# (Arguments otherwise generally the same as studentized version)
gmmq.BS2.test <- function(T.fn,tau,Y,X,Z,Z.excl,Lambda,Lambda.derivative,beta.hat,beta0, Itilde=NULL, Itilde.deriv=NULL, h, structure=c('iid','ts','cluster'),cluster.X.col=NULL,VERBOSE=FALSE,h.adj=1,BREP=99,BLAG=NA) {
  structure <- match.arg(structure)
  if (is.null(Itilde)) Itilde <- get("Itilde.KS17", envir = parent.frame())
  if (is.null(Itilde.deriv)) Itilde.deriv <- get("Itilde.deriv.KS17", envir = parent.frame())
  if (any(is.na(Z))) stop("argument Z should not have any NA entries. Use ret <- gmmq(...,RETURN.Z=TRUE); gmmq.wald.test(...,Z=ret$Z,...)")
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
  # 
  T.hat <- T.fn(beta.hat);  T0 <- T.fn(beta0)
  T.hat.stars <- rep(NA,BREP)
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
    } else stop(sprintf("Argument 'structure' to function gmmq.BS.test() must be 'iid' or 'ts' or 'cluster'"))
    Ystar <- matrix(Y[indstars,],ncol=ncol(Y))
    Xstar <- matrix(X[indstars,],ncol=ncol(X))
    Zexclstar <- matrix(Z.excl[indstars,],ncol=ncol(Z.excl))
    # 
    retstar <- suppressWarnings(tryCatch(gmmq(tau=tau,Y=Ystar,X=Xstar,Z.excl=Zexclstar,dB=length(beta.hat),
                                              Lambda=Lambda, Lambda.derivative=Lambda.derivative,
                                              h=h, VERBOSE=VERBOSE, RETURN.Z=FALSE,
                                              b.init=beta.hat),
                                         error=function(w)NA) )
    if (!is.na(retstar)[1]) T.hat.stars[b] <- T.fn(retstar$b)
  }
  pval.BS <- mean( abs(T.hat.stars-T.hat) > abs(T.hat-T0), na.rm=TRUE)
  if (sum(is.na(T.hat.stars))>0) warning(sprintf("Some NA replications; only %d reps being used. (Increase BREP to increase #non-NA reps.)",sum(!is.na(T.hat.stars))))
  return(list(pval=pval.BS, T.hat=T.hat, BS.T.hats=T.hat.stars))
}





# Estimate asymptotic covariance matrix
# h.adj: should be strictly between 0 and 1/2 if "h=0" (smallest possible bandwidth) was used for estimation, to adjust h to better estimate G; if h.adj=1, then same h is used for G and Sigma estimation (i.e., the value passed in argument h).
cov.est.fn <- function(tau,Y,X,Z,Lambda,Lambda.derivative,beta.hat,Itilde,Itilde.deriv,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE,h.adj=1) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  G.hat <- G.est.fn(Y=Y,X=X,Z=Z,Lambda=Lambda,Lambda.derivative=Lambda.derivative,beta.hat=beta.hat,Itilde.deriv=Itilde.deriv,h=h^h.adj,VERBOSE=VERBOSE)
  Ginv <- tryCatch(solve(G.hat), error=function(w)NA)
  if (is.na(Ginv[1])) return(NA) else {
    LRV.hat <- LRV.est.fn(tau=tau,Y=Y,X=X,Z=Z,Lambda=Lambda,beta.hat=beta.hat,Itilde=Itilde,h=h,structure=structure,cluster.X.col=cluster.X.col,LRV.kernel=LRV.kernel,LRV.ST=LRV.ST,VERBOSE=VERBOSE)
    return(Ginv %*% LRV.hat %*% t(Ginv))
  }
}





# Estimate matrix G
# Z: combines X and Z.excl
# h: not necessarily same as above
G.est.fn <- function(Y,X,Z,Lambda,Lambda.derivative,beta.hat,Itilde.deriv,h,VERBOSE=FALSE) {
  n <- dim(Z)[1]
  L <- Lambda(Y,X,beta.hat)
  Ld <- Lambda.derivative(Y,X,beta.hat)
  tmpsum2 <- t(array(data=Itilde.deriv(-L/h),dim=dim(Z)) * Z) %*% Ld
  return(-tmpsum2/(n*h))
}



# Estimate LRV (matrix Sigma in paper)
# structure: 'iid' for iid data, 'ts' for time series (allowing general weak dependence), 'cluster' for clustered data with independence among clusters but arbitrary dependence among observations within a cluster.
# for structure='cluster': cluster.X.col is a number (like 4) indicating which column of X contains the clustering variable.  (You can always add such a column to X even if Lambda does not use it.)
# LRV.kernel: if structure='ts', should be either `QS' (recommended by Andrews, 1991), 'Bartlett' (like in Newey and West, 1987) or 'uniform' (a.k.a. Truncated in Andrews, 1991); for M-dependent data, can use LRV.kernel='uniform' and LRV.ST=M+1 (e.g. LRV.ST=2 if there is no dependence past 1 lag).
# LRV.ST: if structure='ts', the smoothing parameter S_T from Andrews (1991). For Newey and West (1987), etc., LRV.ST=m+1 where m is how many lags to include in the LRV estimator. Use NA for automated selection.
#
# Kernel functions from Andrews (1991) eqn (2.7)
QS.fn <- function(x) ifelse(x==0,1,(25/(12*pi^2*x^2))*(sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5)))
Bartlett.fn <- function(x) ifelse(abs(x)>=1,0,1-abs(x))
uniform.fn <- function(x) ifelse(abs(x)>=1,0,1) #a.k.a. "Truncated" (Andrews 1991)
# Optimal bandwidths from (5.2) (or 6.2) in Andrews (1991)
QS.ST.fn <- function(alpha2,n) 1.3221*(alpha2*n)^(1/5)
Bartlett.ST.fn <- function(alpha1,n) 1.1447*(alpha1*n)^(1/3)
uniform.ST.fn <- function(alpha2,n) 0.6611*(alpha2*n)^(1/5) #footnote 5
# (6.4) in Andrews (1991)
alpha2.fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/(1-rhos)^8) / sum(ws*sigmas^4/(1-rhos)^4)
alpha1.fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/((1-rhos)^6*(1+rhos)^2)) / sum(ws*sigmas^4/(1-rhos)^4)
# 
LRV.est.fn <- function(tau,Y,X,Z,Lambda,beta.hat,Itilde,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE) {
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  n <- dim(Z)[1]
  if (structure %in% c('iid','ts')) {
    if (structure=='iid') {
      LRV.kernel <- 'uniform'; LRV.lag <- 0; LRV.ST <- 1
    }
    if (missing(LRV.kernel) || !is.character(LRV.kernel)) stop("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS' when structure is 'ts'")
    if (LRV.kernel=='uniform') weight.fn <- uniform.fn else if (LRV.kernel=='Bartlett') weight.fn <- Bartlett.fn else if (LRV.kernel=='QS') weight.fn <- QS.fn else stop(sprintf("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS'; not %s",LRV.kernel))
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
        LRV.ST <- uniform.ST.fn(alpha2=alpha2.fn(rhos=rho.hats,sigmas=sigma.hats),n=n)
      } else if (LRV.kernel=='Bartlett') {
        LRV.ST <- Bartlett.ST.fn(alpha1=alpha1.fn(rhos=rho.hats,sigmas=sigma.hats),n)
      } else if (LRV.kernel=='QS') {
        LRV.ST <- QS.ST.fn(alpha2=alpha2.fn(rhos=rho.hats,sigmas=sigma.hats),n)
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










# EXAMPLE: Job Training Partnership Act data
# gmmq.example.fn() to run
gmmq.example.fn <- function() {
  # Data from Kaplan & Sun (2017) replication .zip: https://drive.google.com/file/d/1N4WmGq6MOxeP5klN3D3iUgGWGj_KH8vP/view
  # or, .csv only: https://drive.google.com/file/d/1AoV-9yqkkzINChmiTorbPk1MV0Du2NxZ/view
  googleID <- "1AoV-9yqkkzINChmiTorbPk1MV0Du2NxZ" # google file ID
  success <-
    tryCatch({jtpa <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", googleID));
    TRUE},  error=function(w) FALSE)
  if (!success) {
    success <-
      tryCatch({df <- read.csv("JTPA_merged.csv"); TRUE},
               error=function(w) FALSE)
    if (!success) {
      stop("Failed to load JTPA_merged.csv from web or local file; should be available at https://drive.google.com/drive/folders/0B-_LUSJVBv20bTNmX3NYeFJKZ2c or through https://kaplandm.github.io")
    }
  }
  ym <- jtpa[jtpa$male==1,c("y")] 
  zm <- jtpa[jtpa$male==1,c("z")]
  dm <- jtpa[jtpa$male==1,c("d")]
  ym <- matrix(ym,5102,1); zm <- matrix(zm,5102,1); dm <- matrix(dm,5102,1)
  om <- matrix(1,nrow(ym),1)
  xm <- jtpa[jtpa$male==1,c(6,7,8,9,10,17,18,12,13,14,15,16,19)]
  xm <- as.matrix(xm)
  LOWh <- 400; HIGHh <- 5e6
  Y <- ym; X <- cbind(om,dm,xm); Z <- cbind(om,zm,xm)
  # 
  tau <- 0.5
  time1 <- system.time(ret1 <- gmmq(tau=tau,Y=cbind(Y,X[,2]),X=X[,-2],Z.excl=matrix(Z[,2],ncol=1),dB=dim(X)[2],Lambda=function(y,x,b)y[,1]-y[,2]*b[1]-x%*%b[-1],Lambda.derivative=function(y,x,b)-cbind(y[,2],x),h=LOWh,VERBOSE=TRUE,RETURN.Z=FALSE,b.init=0))
  time2 <- system.time(ret2 <- gmmq(tau=tau,Y=cbind(Y,X[,2]),X=X[,-2],Z.excl=matrix(Z[,2],ncol=1),dB=dim(X)[2],Lambda=function(y,x,b)y[,1]-y[,2]*b[1]-x%*%b[-1],Lambda.derivative=NULL,h=LOWh,VERBOSE=TRUE))
  print(time1) # w/ derivative = faster
  print(time2) # w/o derivative = slower
  cbind(ret1$b,ret2$b) # identical up to 4 decimals
}
