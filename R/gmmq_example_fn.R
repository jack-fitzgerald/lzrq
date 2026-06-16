# Sourced from David Kaplan's GitHub; see files https://github.com/kaplandm/R/blob/main/ivqr_see.R and https://github.com/kaplandm/R/blob/main/gmmq.R.

gmmq_example_fn <- function() {
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

