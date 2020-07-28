#' @import matlab

kalman_filter <- function(y, A, C, Q, R, Z0, V0) {
  
  
  ###   This function applies a Kalman filter and fixed-interval smoother. The
  ###   script uses the following model:
  ###           y_t = C_t * Z_t + e_t, for e_t ~ N(0, R)
  ###           Z_t = A * Z_{t-1} + mu_t, for mu_t ~ N(0, Q)
  ###   It then applies a fixed-interval smoother using the results from the KF.
  ###  Throughout this file,
  ###    'm' denotes the number of elements in the state vector Z_t.
  ###    'k' denotes the number of elements (observed variables) in y_t.
  ###    'nobs' denotes the number of time periods for which data are observed.
  ###
  ###  Inputs:
  ###    y: k-by-nobs matrix of input data
  ###    A: m-by-m transition matrix 
  ###    C: k-by-m measurement matrix
  ###    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
  ###    R: k-by-k covariance for measurement matrix residuals (e_t)
  ###    Z0: 1-by-m vector, initial value of state
  ###    V0: m-by-m matrix, initial value of factor covariance matrix
  ###
  ###  Outputs:
  ###    Zsmooth: k-by-(nobs+1) matrix, smoothed factor estimates (i.e. Zsmooth[, t + 1] = Z_t|T)
  ###    Vsmooth: k-by-k-by-(nobs+1) array, smoothed factor covariance matrices (i.e. Vsmooth[, , t + 1) = Cov(Z_t|T))
  ###    VVsmooth: k-by-k-by-nobs array, lag 1 factor covariance matrices (i.e. Cov(Z_t, Z_t-1|T))
  ###    loglik: scalar, log-likelihood
  ###
  ### This version: Fernando Cantu, 2020-06-01
  ###
  
  ###
  ### Apply Kalman filter
  ###
  
  ### Initialise
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
    output <- remove_missing(as.matrix(array(y[, t])), C, R)
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
  
  return(list(Zsmooth = ZmT, Vsmooth = VmT, VVsmooth = VmT_lag, loglik = loglik))
  
}   