lzrq = function(formula, data, tau = 0.5, floor_val = -1e300, ...) {
  
  ##################
  ##### ERRORS #####
  ##################
  
  #If quantreg is not installed...
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    
    #... then stop the function
    stop("package 'quantreg' is required; install with: install.packages('quantreg')")
    
  }
  
  ###################
  ##### SETUP #######
  ###################
  
  #Extract model frame (handles NA removal, data lookup, subset, weights, etc.)
  mf = model.frame(formula, data = data, na.action = na.omit)
  
  #Extract raw outcome (left-hand side, before log transformation)
  Y = model.response(mf)
  
  #If Y contains no strictly positive values...
  if (!any(Y > 0)) {
    
    #... then stop the function
    stop("'Y' has no positive values; log quantile regression requires at least some positive outcomes")
    
  }
  
  #Number of observations and non-positive values
  n        = length(Y)
  n_nonpos = sum(Y <= 0)
  
  #Smallest positive value of the outcome
  ln_min_pos = log(min(Y[Y > 0]))
  
  #Step 1: initial sentinel = midpoint between floor and ln(min positive Y)
  sentinel = (floor_val + ln_min_pos) / 2
  
  ###########################
  ##### ITERATIVE BISECT ####
  ###########################
  
  success = FALSE
  
  for (iter in 1:10) {
    
    #Build transformed outcome: ln(Y) for Y > 0, sentinel for Y <= 0
    logY = ifelse(Y > 0, log(Y), sentinel)
    
    #Update the model frame with the transformed outcome
    mf[[1]] = logY
    
    #Run quantile regression on transformed outcome, passing all extra args through
    fit = quantreg::rq(formula = formula(mf), data = mf, tau = tau, ...)
    
    #Check whether all fitted values exceed the sentinel
    yhat = fit$fitted.values
    
    if (all(yhat > sentinel)) {
      
      success = TRUE
      break
      
    }
    
    #Not all fitted values above sentinel: move sentinel halfway to floor
    sentinel = (floor_val + sentinel) / 2
    
  }
  
  #####################
  ##### CONVERGENCE ###
  #####################
  
  #If convergence failed after 10 iterations...
  if (!success) {
    
    #... then stop the function
    stop(paste0(
      "Convergence failure: after 10 iterations, some fitted values from quantile regression on the transformed outcome remain at or below the psi value; results suppressed"
    ))
    
  }
  
  ###################
  ##### OUTPUT ######
  ###################
  
  #Attach lzrq-specific fields to the rq object
  fit$sentinel = sentinel
  fit$n_nonpos = n_nonpos
  
  #Store original call separately for display only; leave fit$call as internal rq call
  #so that all rq methods (predict, etc.) keep working
  fit$call_lzrq      = match.call()
  fit$call_lzrq[[1]] = as.name("lzrq")
  
  #Set class to inherit from rq so all rq methods work automatically
  class(fit) = c("lzrq", "rq")
  
  #Return output invisibly (print.lzrq handles display)
  return(fit)
  
}

#####################
##### S3 METHODS ####
#####################

print.lzrq = function(x, ...) {
  cat(sprintf("\nQuantile regression (log outcome)     Number of obs        = %8d\n", length(x$fitted.values)))
  cat(sprintf("Outcome: log(%s)      Number non-positive  = %d\n\n", deparse(x$call_lzrq$formula[[2]]), x$n_nonpos))
  #Temporarily swap in the lzrq call for display, then restore
  rq_call    = x$call
  x$call     = x$call_lzrq
  NextMethod(x, ...)
  message("\nPlease cite the papers underlying this command:")
  message("  Fitzgerald, J., Adema, J., Fiala, L., Kujansuu, E., & Valenta, D. (2026). Non-Robustness in Log-Like Specifications. MetaArXiv. https://doi.org/10.31222/osf.io/juda7_v1")
  message("  Liu, X., & Kaplan, D. M. (2025). Quantile Regression with Log(0) Outcomes. https://drive.google.com/file/d/1F3dnhm8MrlO5aRrGt48rBWAEaBqdCBH-/view")
  invisible(x)
}

summary.lzrq = function(object, se = "nid", covariance = FALSE, wald = FALSE, ...) {
  cat(sprintf("\nQuantile regression (log outcome)      Number of obs        = %8d\n", length(object$fitted.values)))
  cat(sprintf("Outcome: log(%s)      Number non-positive  = %d\n\n", deparse(object$call_lzrq$formula[[2]]), object$n_nonpos))
  rq_call     = object$call
  object$call = object$call_lzrq
  result = NextMethod()
  print(result)
  message("\nPlease cite the papers underlying this command:")
  message("  Fitzgerald, J., Adema, J., Fiala, L., Kujansuu, E., & Valenta, D. (2026). Non-Robustness in Log-Like Specifications. MetaArXiv. https://doi.org/10.31222/osf.io/juda7_v1")
  message("  Liu, X., & Kaplan, D. M. (2025). Quantile Regression with Log(0) Outcomes. https://drive.google.com/file/d/1F3dnhm8MrlO5aRrGt48rBWAEaBqdCBH-/view")
  invisible(result)
}