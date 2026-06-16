# Implementation of estimator from
# Kaplan and Sun (2017), "Smoothed estimating equations for instrumental variables quantile regression" (Econometric Theory)
# Also available in Stata; visit  https://kaplandm.github.io/#sivqr
# Questions? Comments? email Dave  (https://kaplandm.github.io)
# 
# ! an intercept will be included automatically;
#   you don't need one in X.exog or Z.excl
# ! careful of the order of coefficients: first endogenous,
#   then exogenous, then intercept
# ! standard errors assume iid sampling
# return value "h" is the actual bandwidth used;
#    "hhat" is the plug-in (or desired) bandwidth
# As you can see, ivqr.see() is basically a wrapper for gmmq()

# Load code
prob.loaded <- exists("gmmq")
success <-
  tryCatch({source('https://raw.githubusercontent.com/kaplandm/R/main/gmmq.R');TRUE},
           error=function(w) FALSE)
if (!success) {
  success <- tryCatch({source("gmmq.R");TRUE}, error=function(w) FALSE)
  if (!success) {
    if (prob.loaded) {
      warning("Couldn't load gmmq.R, but it seems like you already did.")
    } else {
      stop("Failed to source() gmmq.R from web or local file (possibly because a package like pracma is not installed? check the other warnings/errors).  You may download and source() it yourself, or at least make sure it's in your getwd().  Currently available at https://github.com/kaplandm/R")
    }
  }
}

# tau: 0<tau<1 is the quantile index, like tau=0.5 for median
# Y: a vector of the observed scalar outcome variable for each observation
# D: a matrix of the endogenous regressors; each row is a different observation, 
#    each column is a different variable
# X.exog: a matrix similar to D but with exogenous regressors
# Z.excl: a matrix like D or X.exog but with the excluded instruments (exogenous 
#         variables that are not part of the model for Y); the variables in 
#         Z.excl instrument for the endogenous regressors in D
# h: desired smoothing bandwidth, or else plug-in bandwidth is used if h is not
#    specified
# b.init: initial value of coefficient vector for numerical search, or else
#         QR is used for initial value if b.init is not specified
# iidSE: report standard errors? (TRUE/FALSE) only valid for iid sampling
ivqr.see <- function(tau, Y, D=NULL, X.exog=NULL, Z.excl=NULL, h, b.init=NULL, iidSE=TRUE) {
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
    hhat <- ivqr.bw(p=tau, Y=Y, X=cbind(D,X), 
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

# Implementation of plug-in bandwidth from
# Kaplan and Sun (2017), "Smoothed estimating equations for instrumental variables quantile regression" (Econometric Theory)
# Questions? Comments? email Dave: kaplandm@missouri.edu
# 
# Use gmmq() with h=0 to get b.init
# Then use ivqr.bw() to compute the "optimal" bandwidth
# Then re-run gmmq() with that h
# 
# p: quantile index, 0<p<1 (0.5=median, etc.)
# Y: vector of scalar outcome observations
# X: matrix of regressors
# b.init: initial estimate
# So, residuals are Y-X%*%b.init

ivqr.bw <- function(p,Y,X,b.init) {
  r <- 4
  if (!is.matrix(Y)) Y <- matrix(Y,ncol=1)
  # Defaults, initialization
  n <- length(Y) #number of observations, a.k.a. sample size
  d <- ncol(X) #number of regressors
  
  # 3-parameter student's t PDF
  dt.ls.fn <- function(x,m,sigma,df) dt(x=(x-m)/sigma,df=df)/sigma
  # GEV density function (PDF)
  dgev.ls.fn <- function(x,xi,mu=0,sig=1) {
    if (sig<=0) {
      return(rep(0,length(x)))
    } else if (xi!=0) {
      tmp <- (1 + (xi * (x - mu))/sig)
      return(ifelse(tmp>0,(tmp^(-1/xi - 1)*exp(-tmp^(-1/xi)))/sig,0))
    } else {
      tmp <- exp(-(x-mu)/sig)
      return((tmp^(xi+1)*exp(-tmp))/sig)
    }
  }
  
  # For different possible parametric distributions, 
  # generate the objects needed for the bandwidth as
  # functions of distribution parameters.
  fN0_fn <- function(nmu,sig) {
    dnorm(0,nmu,sig)
  }
  fN0r1_fn <- function(nmu,sig) {
    dnorm(0,nmu,sig)%*%(nmu/sig^4)%*%(-3+(nmu%*%nmu)/sig^2)
  } 
  ft0_fn <- function(v,tmu,tsig) {
    if (v>100) {
      dnorm(0,tmu,tsig)
    } else {
      pi^(-1/2)%*%((1/sqrt(v)*gamma(v*(1/2)+1/2)*((tmu^2*1/tsig^2)/v+1)^(v*(-1/2)-1/2))/(tsig*gamma(v*(1/2))))
    }
  }
  ft0r1_fn <- function(v,tmu,tsig) {
    if (v>100) {
      fN0r1_fn(tmu,tsig)
    } else {
      pi^(-1/2)%*%((tmu*1/tsig^5*1/v^(5/2)*gamma(v*(1/2)+1/2)*(v*(1/2)+1/2)*(v*(1/2)+3/2)*((tmu^2*1/tsig^2)/v+1)^(v*(-1/2)-5/2)*-12)/gamma(v*(1/2))+(tmu^3*1/tsig^7*1/v^(7/2)*gamma(v*(1/2)+1/2)*(v*(1/2)+1/2)*(v*(1/2)+3/2)*(v*(1/2)+5/2)*((tmu^2*1/tsig^2)/v+1)^(v*(-1/2)-7/2)*8)/gamma(v*(1/2)))
    }
  }
  fgam0_fn <- function(k,theta,Ueval) {
    dgamma(Ueval,k,scale=theta)
  }
  fgam0r1_fn <- function(k,theta,Ueval) {
    # ORIGINAL (equivalent but get Inf "errors" b/c not grouping Ueval/theta): -(Ueval^(k-1)*1/theta^3*theta^(-k)*exp(-Ueval/theta))/gamma(k)+(Ueval^(k-2)*1/theta^2*theta^(-k)*exp(-Ueval/theta)*(k-1)*3)/gamma(k)-(Ueval^(k-3)*theta^(-k)*exp(-Ueval/theta)*(k-1)*(k-2)*3)/(theta*gamma(k))+(Ueval^(k-4)*theta^(-k)*exp(-Ueval/theta)*(k-1)*(k-2)*(k-3))/gamma(k)
    -((Ueval/theta)^(k-1)*1/theta^3*theta^(-1)*exp(-Ueval/theta))/gamma(k)+((Ueval/theta)^(k-2)*1/theta^2*theta^(-2)*exp(-Ueval/theta)*(k-1)*3)/gamma(k)-((Ueval/theta)^(k-3)*theta^(-3)*exp(-Ueval/theta)*(k-1)*(k-2)*3)/(theta*gamma(k))+((Ueval/theta)^(k-4)*theta^(-4)*exp(-Ueval/theta)*(k-1)*(k-2)*(k-3))/gamma(k)
  }
  fGEV0_fn <- function(xi,sig,gmu) {
    dgev.ls.fn(0,xi,gmu,sig)
  } 
  fGEV0r1_fn <- function(xi,sig,gmu) {
    1/sig^4*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(-(gmu*xi)/sig+1)^(-3/xi-3)-1/sig^4*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^xi*(-(gmu*xi)/sig+1)^(-3/xi-3)*(xi+1)-1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-2/xi-3)*2-1/sig^4*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^xi*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-2/xi-2)*(xi+1)*2+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-1)*(-(gmu*xi)/sig+1)^(-3/xi-3)*(xi+1)+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^xi*(-(gmu*xi)/sig+1)^(-2/xi-3)*(2/xi+2)*(xi+1)*2+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-1)*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-2/xi-2)*(xi+1)*2-1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-1)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-2/xi-3)*(xi+1)*2-1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-1/xi-2)-1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-2)*(-(gmu*xi)/sig+1)^(-3/xi-3)*(xi-1)*(xi+1)+1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(1/xi+1)*(1/xi+2)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-3.0)+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1.0)^(-1.0/xi))*((-(gmu*xi)/sig+1.0)^(-1.0/xi))^xi*(1.0/xi+1.0)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-1.0)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-2.0)*(xi+1.0)*2.0-1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1.0)^(-1.0/xi))*((-(gmu*xi)/sig+1.0)^(-1.0/xi))^xi*(1.0/xi+1.0)*(1.0/xi+2.0)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-3.0)*(xi+1.0)-1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1.0)^(-1.0/xi))*((-(gmu*xi)/sig+1.0)^(-1.0/xi))^(xi-1.0)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-1/xi-2)*(xi+1)
  }
  hopt_estfn1_t <- function(Gsqiv,CKv,mu,sig,v,d,n) {
    ((factorial(r)^2*(1-Gsqiv)*ft0_fn(v,mu,sig)*d)/(2*r*CKv^2*ft0r1_fn(v,mu,sig)^2*n))^(1/(2*r-1))
  }
  hopt_estfn1_N <- function(Gsqiv,CKv,mu,sig,d,n) {
    ((factorial(r)^2*(1-Gsqiv)*fN0_fn(mu,sig)*d)/(2*r*CKv^2*fN0r1_fn(mu,sig)^2*n))^(1/(2*r-1))
  }
  hopt_estfn1_gam <- function(Gsqiv,CKv,k,theta,Ueval,d,n) {
    ((factorial(r)^2*(1-Gsqiv)*fgam0_fn(k,theta,Ueval)*d)/(2*r*CKv^2*fgam0r1_fn(k,theta,Ueval)^2*n))^(1/(2*r-1))
  }
  hopt_estfn1_GEV <- function(Gsqiv,CKv,xi,sig,mu,d,n) {
    ((factorial(r)^2*(1-Gsqiv)*fGEV0_fn(xi,sig,mu)*d)/(2*r*CKv^2*fGEV0r1_fn(xi,sig,mu)^2*n))^(1/(2*r-1))
  }
  # # Set smoothing fn G() and its derviative G'()
  # G_symfn <- function(u) {
  #   ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0))
  # } 
  # Gp_symfn <- function(u) {
  #   ifelse(u > -1 & u < 1, (105/64)*(1-5*u^2+7*u^4-3*u^6), 0)
  # } 
  Gsqint_val <- 394/429
  CK_val <- -1/33
  
  # # Add bandwidth as second argument to G(), G'()
  # Gfn <- function(v,h) { G_symfn(v/h) }
  # Gpfn <- function(v,h) { Gp_symfn(v/h) }
  
  # Plug-in bandwidth as function of distributional parameters, for different possible parametric distributions.
  hopt_estfn_t <- function(mu,sig,v,d,n) {
    hopt_estfn1_t(Gsqint_val,CK_val,mu,sig,v,d,n)
  }
  hopt_estfn_N <- function(mu,sig,d,n) {
    hopt_estfn1_N(Gsqint_val,CK_val,mu,sig,d,n)
  }
  hopt_estfn_gam <- function(k,theta,p,d,n) {
    hopt_estfn1_gam(Gsqint_val,CK_val,k,theta,qgamma(p,k,scale=theta),d,n)
  }
  hopt_estfn_GEV <- function(xi,sig,mu,d,n) {
    hopt_estfn1_GEV(Gsqint_val,CK_val,xi,sig,mu,d,n)
  }
  
  Uhat <- Y - X%*%b.init
  # Shift Uhat dist'n such that p-quantile=0
  Uhat <- Uhat - quantile(Uhat, p, na.rm=TRUE)
  if (!require(MASS)) {
    stop('Please install the R package MASS to run this code.')
  }
  t.MLE.est <- tryCatch(fitdistr(x=Uhat,densfun=dt.ls.fn,start=list(m = mean(Uhat), sigma = sd(Uhat), df = 100),lower=c(-Inf,0,1), control=list(maxit=200)),
                        warning=function(w) NA,
                        error=function(w) NA
  ) # If warning/error, don't use (set to Inf)
  if (is.na(t.MLE.est[1])) {
    tnlogL <- Inf
    hoptests_t <- Inf
  } else {
    tnlogL <- -t.MLE.est$loglik
    v <- max(1,round(t.MLE.est$estimate["df"]))
    mu <- t.MLE.est$estimate["m"]
    sig <- t.MLE.est$estimate["sigma"]
    hoptests_t <- hopt_estfn_t((mu+(0==mu)*.01),sig,v,d,n)
  }
  tmp.sig <- sd(Uhat)*sqrt(6)/pi
  tmp.mu <- mean(Uhat) - tmp.sig*(-digamma(1)) #approx 0.5772... (Euler's Constant)
  gev.MLE.est <- tryCatch(suppressWarnings(fitdistr(x=Uhat,densfun=dgev.ls.fn,list(xi=0, mu=tmp.mu, sig=2*tmp.sig), method='BFGS', control=list(maxit=200))), 
                          #warning=function(w) NA,
                          error=function(w) {
                            tryCatch(suppressWarnings(fitdistr(x=Uhat,densfun=dgev.ls.fn, list(xi=0, mu=tmp.mu, sig=tmp.sig), lower=c(-Inf,-Inf,0), control=list(maxit=200))), 
                                     error=function(w2)NA)
                          }
  )
  if(is.na(gev.MLE.est[1]) || any(!is.finite(gev.MLE.est$estimate))) {
    gevnlogL <- Inf
    hoptests_gev <- Inf
  } else {
    gevnlogL <- -gev.MLE.est$loglik
    xi <- gev.MLE.est$estimate["xi"]
    sig <- gev.MLE.est$estimate["sig"]
    mu <- gev.MLE.est$estimate["mu"]
    hoptests_gev <- hopt_estfn_GEV(xi=xi,sig=sig,mu=mu,d=d,n=n)
  }
  gam.MLE.est <- tryCatch(suppressWarnings(fitdistr(x=Uhat-min(Uhat)+.Machine$double.eps,densfun=dgamma,start=list(shape=1,scale=sd(Uhat)),lower=.Machine$double.eps, control=list(maxit=200))),
                          # warning=function(w) NA,
                          error=function(w) NA
  )
  if(is.na(gam.MLE.est[1]) || any(!is.finite(gam.MLE.est$estimate))) {
    gamnlogL <- Inf
    hoptests_gam <- Inf
  } else {
    gamnlogL <- -gam.MLE.est$loglik
    k <- gam.MLE.est$estimate["shape"]
    theta <- gam.MLE.est$estimate["scale"]
    hoptests_gam <- hopt_estfn_gam(k,theta,p,d,n)
  }
  mu <- mean(Uhat); sig <- sd(Uhat)
  hoptests_N <- hopt_estfn_N((mu+(0==mu)*.01),sig,d,n)
  h <- min(hoptests_t,hoptests_N,hoptests_gev,hoptests_gam , na.rm=TRUE)
  return(c(h=h))
}