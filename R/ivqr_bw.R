# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

ivqr_bw <- function(p,Y,X,b.init) {
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
