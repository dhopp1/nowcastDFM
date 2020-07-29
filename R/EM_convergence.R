#' @importFrom pracma eps

EM_convergence <- function(loglik, prev_loglik, threshold) {
  
  ###
  ###  This function checks whether EM has converged. Convergence occurs if the slope of the
  ###  log-likelihood function falls below 'threshold' (i.e. f(t) - f(t-1)| / avg < threshold) 
  ###  where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log likelihood at iteration t. 
  ###  This stopping criterion is from Numerical Recipes in C (pg. 423). With MAP estimation (using priors),
  ###  the likelihood can decrease even if the mode of the posterior increases.
  ###
  ###  Inputs:
  ###     - loglik:        log-likelihood from current EM iteration
  ###     - prev_loglik:   log-likelihood from previous EM iteration
  ###     - threshold:     convergence threshhold
  ###
  ###  Outputs:
  ###     - converged:   1 if convergence criteria satisfied, 0 otherwise
  ###     - decrease:    1 if loglikelihood decreased, 0 otherwise
  ###
  ###  This version: Fernando Cantu, 2020-07-03
  ###
  
  ### Initialize output
  converged <- 0
  decrease <- 0
  
  ### Check if log-likelihood decreases
  if ((loglik - prev_loglik) < -1e-3) {    # Allows for a little imprecision
    decrease <- 1
  }
  
  ### Check convergence criteria
  delta <- abs(loglik - prev_loglik)
  avg_loglik = (abs(loglik) + abs(prev_loglik) + eps()) / 2
  if ((delta / avg_loglik) < threshold){
    converged <- 1
  }
  
  return(list(converged = converged, decrease = decrease))
  
}