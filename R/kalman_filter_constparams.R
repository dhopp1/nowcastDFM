#' @importFrom pracma inv pinv

kalman_filter_constparams <- function(data, params, lag) {
  
  
  ###   This function applies a Kalman filter for news calculation step, when model parameters
  ###   are already estimated. This procedure only smoothes and fills missing data for a given data matrix
  ###
  ###  Inputs:
  ###    data:  input data matrix (transformed, but non-standardized)
  ###    params:  parameters, results structure from previously estimated DFM
  ###    lag:  number of lags
  ###
  ###  Outputs:
  ###    Plag:  Smoothed factor covariance for transition matrix
  ###    Vsmooth:  Smoothed factor covariance matrix
  ###    X_smooth: Smoothed data matrix
  ###    F:  Smoothed factors
  ###
  ### This version: Fernando Cantu, 2020-07-07
  ###
  
  ###
  ### Apply Kalman filter
  ###
  
  ### Initialise
  Z0 <- params$Z0
  V0 <- params$V0
  A <- params$A
  C <- params$C
  Q <- params$Q
  R <- params$R
  means <- params$means
  sdevs <- params$sdevs
  T <- dim(data)[1]
  y <- t((data - repmat(as.numeric(means), T, 1)) / repmat(as.numeric(sdevs), T, 1))
  
  m <- dim(C)[2]
  nobs <- dim(y)[2]
  Zm  <- array(NA, dim = c(m, nobs))         # Z_t | t-1 (prior)
  Vm  <- array(NA, dim= c(m, m, nobs))       # V_t | t-1 (prior)
  ZmU <- array(NA, dim = c(m, nobs + 1))     # Z_t | t (posterior/updated)
  VmU <- array(NA, dim = c(m, m, nobs+1))    # V_t | t (posterior/updated)
  ZmT <- array(0, dim = c(m, nobs+1))        # Z_t | T (smoothed states)
  VmT <- array(0, dim = c(m, m, nobs + 1))   # V_t | T = Cov(Z_t|T) (smoothed factor covariance)
  VmT_lag <- array(NA, dim = c(m, m, nobs))   # Cov(Z_t, Z_t-1|T) (smoothed lag 1 factor covariance)
  loglik <- 0
  
  ### Initial values
  
  Zu <- Z0    # Z_0|0 (In loop, Zu gives Z_t | t)
  Vu <- V0    # V_0|0 (In loop, Vu gives V_t | t)
  ZmU[, 1] <- Zu
  VmU[, , 1] <- Vu
  
  ### Kalman filter
  for (t in 1:nobs) {
    
    #-# Calculate prior distribution
    # Use transition eqn to create prior estimate for factor, i.e. Z = Z_t|t-1
    Z <- A %*% Zu
    # Prior covariance matrix of Z (i.e. V = V_t|t-1):
    # Var(Z) = Var(A*Z + u_t) = Var(A*Z) + Var(\epsilon) = A*Vu*A' + Q
    V <- A %*% Vu %*% t(A) + Q
    V <- 0.5 * (V + t(V))  # Trick to make symmetric
    
    #-# Calculate posterior disttribution
    # Remove missing series: These are removed from Y, C, and R
    output <- remove_missing(y[, t], C, R)
    yt <- output$yt; Ct <- output$C; Rt <- output$R
    if (length(Ct) == m) Ct <- t(Ct)
    # Check if yt contains no data; if this is the case, replace Zu and Vu with prior
    if (length(yt) == 0) {
      Zu <- Z
      Vu <- V
    } else {
      # Steps for variance and population regression coefficients:
      # Var(c_t * Z_t + e_t) = c_t * Var(A) * c_t' + Var(u) = c_t * V * c_t' + R
      VC <- V %*% t(Ct)  
      iF  <- inv(Ct %*% VC + Rt)
      # Matrix of population regression coefficients (QuantEcon eqn #4)
      VCF <- VC %*% iF
      # Gives difference between actual and predicted measurement matrix values
      innov <- yt - Ct %*% Z
      # Update estimate of factor values (posterior)
      Zu  = Z  + VCF %*% innov;
      # Update covariance matrix (posterior) for time t
      Vu <- V - VCF %*% t(VC)
      Vu <- 0.5 * (Vu + t(Vu))  # Trick to make symmetric
      # Update log likelihood 
      loglik = loglik + 0.5 * (log(det(iF)) - t(innov) %*% iF %*% innov)
    }
    
    ### Store output
    # Store covariance and observation values for t-1 (priors)
    Zm[, t] <- Z
    Vm[, , t] <- V
    # Store covariance and state values for t (posteriors), i.e. Zu = Z_t|t & Vu = V_t|t
    ZmU[, t + 1]  <- Zu
    VmU[, , t + 1] <- Vu
    
  }
  
  ### Store Kalman gain k_t
  if (length(yt) == 0) {
    k_t <- zeros(m, m)
  } else {
    k_t <- VCF %*% Ct
  }
  
  
  ###
  ### Apply fixed interval smoother
  ###
  
  # Fill the final period of ZmT & VmT with posterior values from KF
  ZmT[, nobs + 1] <- drop(ZmU[, nobs + 1])
  VmT[, , nobs + 1] <- drop(VmU[, , nobs + 1])
  # Initialize VmT_1 lag 1 covariance matrix for final period
  VmT_lag[, , nobs] <- (eye(m) - k_t) %*% A %*% drop(VmU[, , nobs])
  # Used for recursion process, see companion file for details
  J_2 <- drop(VmU[, , nobs]) %*% t(A) %*% pinv(drop(Vm[, , nobs]))
  
  ### Run smoothign algorithm
  # Loop through time reverse-chronologically (starting at final period nobs)
  for (t in nobs:1) {
    
    # Store posterior and prior factor covariance values 
    VmUt <- drop(VmU[, , t])
    Vmt <- drop(Vm[, , t])
    # Store previous period smoothed factor covariance and lag-1 covariance
    Vt <- drop(VmT[, , t + 1])
    Vt_lag <- drop(VmT_lag[, , t])
    J_1 <- J_2
    # Update smoothed factor estimate
    ZmT[, t] <- ZmU[, t] + J_1 %*% (ZmT[, t + 1] - A %*% ZmU[, t]) 
    # Update smoothed factor covariance matrix
    VmT[, , t] <- VmUt + J_1 %*% (Vt - Vmt) %*% t(J_1)
    if (t > 1) {
      # Update weight
      J_2 <- drop(VmU[, , t - 1]) %*% t(A) %*% pinv(drop(Vm[, , t - 1]))
      # Update lag 1 factor covariance matrix 
      VmT_lag[, , t - 1] <- VmUt %*% t(J_2) + J_1 %*% (Vt_lag - A %*% VmUt) %*% t(J_2)
    }
    
  }
  
  ###
  ### Prepare output
  ###
  
  Vs <- VmT[, , -1]  # Smoothed factor covariance for transition matrix
  Vf <- VmU[, , -1]  # Filtered factor posterior covariance
  Zsmooth <- ZmT     # Smoothed factors
  Vsmooth <- VmT     # Smoothed covariance value
  Plag <- list()
  Plag[[1]] <- Vs
  
  if (lag > 0) {
    for (jk in 1:lag) {
      Plag[[jk + 1]] <- array(NA, dim = c(m, m, nobs)) 
      for (jt in dim(Plag[[1]])[3]:(lag + 1)) {
        As <- Vf[, , jt - jk] %*% t(A) %*% pinv(A %*% Vf[, , jt - jk] %*% t(A) + Q)
        Plag[[jk + 1]][,  , jt] <- As %*% Plag[[jk]][,  , jt]
      }
    }
  }
  
  Zsmooth <- t(Zsmooth)
  x_sm <- Zsmooth[-1, ] %*% t(C)   # Factors to series representation
  X_smooth <- repmat(as.numeric(sdevs), T, 1) * x_sm + repmat(as.numeric(means), T, 1)   # Standardized to unstandardized
  
  return(list(Plag = Plag, Vsmooth = Vsmooth, X_smooth = X_smooth, F = Zsmooth[-1, ]))
  
}   