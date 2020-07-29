#' @importFrom pracma inv kron blkdiag

init_conds <- function(X, p, blocks, R_mat, q, nM, nQ, index_freq) {
  ###
  ### This function creates initial parameter estimates. These are initial inputs in the EM algorithm, 
  ### which re-estimates these parameters using the Kalman filter
  ###
  ### Inputs:
  ###    - X: Standardised data
  ###    - p: Number of lags in transition equation
  ###    - blocks: Block structure
  ###    - nM: Number of monthly variables
  ###    - nQ: Number of quarterly variables
  ###    - index_freq: Indicator for monthly (1) and quarterly variables (0)
  ###
  ### Outputs:
  ###    - A:   Transition matrix
  ###    - C:   Measurement matrix
  ###    - Q:   Covariance for transition equation residuals
  ###    - R:   Covariance for measurement equation residuals
  ###    - Z0:  Initial value of state
  ###    - V0:  Initial value of covariance matrix
  ###
  ### This version: Fernando Cantu, 2020-05-31
  ###
  
  num_blocks <- ncol(blocks)
  
  output <- fill_na(X)
  data_full <- output$data_full
  ind_na <- output$ind_na
  T <- nrow(data_full)
  N <- ncol(data_full)
  
  data_temp <- data_full
  data_temp[ind_na] <- NA
  res <- data_full
  res_temp <- data_temp
  
  pC <- 5
  ppC <- max(p, pC)
  C = NULL;
  A = NULL;
  Q = NULL;
  V0 = NULL;
  
  ind_na[1:(pC - 1), ] <- TRUE
  
  for (i in 1:num_blocks) {
    
    ### Measurement equation ------------------------------------------------------------
    C_i <- zeros(ncol(data_full), ppC)   # Initialize state variable matrix helper
    idx_i <- which(blocks[, i] == 1)    # Indices of series loading in block i
    idx_iM <- idx_i[idx_i <= nM]  # Monthly series indices that load in block i
    idx_iQ <- idx_i[idx_i > nM]    # Quarterly series indices that loaded in block i
    eigendecomp <- eigen(cov(res[, idx_iM]))   # Eigenvector of cov matrix of monthly data and largest eigenvalue
    d <- eigendecomp$values[1]   # Largest eigenvalue...
    v <- eigendecomp$vectors[, 1]   # ...and the associated eigenvector
    if(sum(v) < 0) v <- -v  # Flip sign for clearer output (not required, but facilitates reading)
    C_i[idx_iM, 1] <- v
    f <- as.matrix(res[, idx_iM]) %*% v   # Data projection for eigenvector direction
    F <- NULL   # Lag matrix
    for (j in 0:(max(p+1, pC) - 1)) {
      F <- cbind(F, f[(pC-j):(nrow(f)-j)])
    }
    ff <- F[, 1:pC]    # Projected data with lag structure, so pC-1 fewer entries
    for (j in idx_iQ) {   #  Loop for quarterly variables
      xx_j <- res_temp[pC:nrow(res_temp), j]   # For series j, values are dropped to accommodate lag structure
      if (sum(!is.na(xx_j)) < size(ff, 2) + 2) {
        xx_j <- res[pC:nrow(res), j]    # Replaces xx_j with spline if too many NaNs
      }
      ff_j <- ff[!is.na(xx_j), ]
      xx_j <- xx_j[!is.na(xx_j)]
      iff_j <- inv(t(ff_j) %*% ff_j)
      Cc <- iff_j %*% t(ff_j) %*% xx_j   # OLS
      Cc <- Cc - iff_j %*% t(R_mat) %*% inv(R_mat %*% iff_j %*% t(R_mat)) %*% (R_mat %*% Cc - q)   # Spline data monthly to quarterly conversion
      C_i[j, 1:pC] <- t(Cc)   # Place in output matrix
    }
    ff <- rbind(zeros(pC-1, pC), ff)   # Zeros in first pC-1 entries (replacing dropped from lag)
    res <- res - ff %*% t(C_i)   # Residual calculations
    res_temp <- res
    res_temp[ind_na] <- NA
    C <- cbind(C, C_i)   # Combine past loadings together
    
    ### Transition equation ------------------------------------------------------------
    z <- F[, 1]   # Projected data (no lag)
    Z <- F[, 2:(p+1)]   # Data with lag 1
    A_i <- zeros(ppC, ppC)   # Initialize transition matrix
    A_temp <- inv(t(Z) %*% Z) %*% t(Z) %*% z    # OLS: gives coefficient value AR(p) process
    A_i[1, 1:p] <- t(A_temp)
    A_i[2:nrow(A_i), 1:(ppC-1)] <- eye(ppC-1)
    ###
    Q_i <- zeros(ppC, ppC)
    e <- z - Z %*% A_temp    # VAR residuals
    Q_i[1, 1] <- cov(e)   # VAR covariance matrix
    initV_i <- matrix(inv(eye(ppC ^ 2) - kron(A_i, A_i)) %*% as.vector(Q_i), nrow = ppC, ncol = ppC)
    # Gives top left block for the transition matrix
    if (is.null(A)) {A <- A_i} else {A <- blkdiag(A, A_i)}
    if (is.null(Q)) {Q <- Q_i} else {Q <- blkdiag(Q, Q_i)}
    if (is.null(V0)) {V0 <- initV_i} else {V0 <- blkdiag(V0, initV_i)}
  }
  
  eyeN <- eye(N)   # Used inside measurement matrix
  eyeN <- eyeN[, ifelse(index_freq == 1, TRUE, FALSE)]
  
  C <- cbind(C, eyeN)
  C <- cbind(C, rbind(zeros(nM, 5 * nQ), t(kron(eye(nQ), c(1, 2, 3, 2, 1)))))  # Monthly-quarterly agreggation scheme
  R <- diag(apply(res_temp, 2, var, na.rm = T))  # Initialize covariance matrix for transition matrix
  
  ii_idio <- which(ifelse(index_freq == 1, TRUE, FALSE))    # Indicies for monthly variables
  n_idio <- length(ii_idio)  # Number of monthly variables
  BM <- zeros(n_idio)        # Initialize monthly transition matrix values
  SM <- zeros(n_idio)        # Initialize monthly residual covariance matrix values
  
  for (i in 1:n_idio) {   # Loop for monthly variables
    # Set measurement equation residual covariance matrix diagonal
    R[ii_idio[i], ii_idio[i]] <- 1e-04
    # Subsetting series residuals for series i
    res_i <- res_temp[, ii_idio[i]]
    # Returns number of leading/ending zeros
    leadZero <- max(which(1:T == cumsum(is.na(res_i))))
    endZero <- which(1:T == cumsum(is.na(rev(res_i))))
    endZero <- ifelse(length(endZero) == 0, 0, max(endZero))
    # Truncate leading and ending zeros
    res_i <- res[, ii_idio[i]]
    res_i <- res_i[(leadZero + 1):(length(res_i) - endZero)]
    # Linear regression: AR(1) process for monthly series residuals
    BM[i, i] <- inv(res_i[1:(length(res_i)-1)] %*% res_i[1:(length(res_i)-1)]) %*% 
      res_i[1:(length(res_i)-1)] %*% res_i[2:length(res_i)]
    SM[i, i] <- var(res_i[2:length(res_i)] - res_i[1:(length(res_i)-1)] * BM[i, i])    # Residual covariance matrix
  }
  
  Rdiag <- diag(R)
  sig_e <- Rdiag[(nM + 1):N] / 19
  Rdiag[(nM+1):N] <- 1e-04
  R <- diag(Rdiag)    # Covariance for obs matrix residuals
  
  # For BQ, SQ
  rho0 <- 0.1
  temp <- zeros(5)
  temp[1, 1] <- 1
  
  # Blocks for covariance matrices
  if (length(sig_e) == 1) {
    SQ <- kron(((1-rho0^2)*sig_e), temp)
  } else {
    SQ <- kron(diag((1-rho0^2)*sig_e), temp)
  }
  BQ <- kron(eye(nQ), rbind(c(rho0, zeros(1, 4)), cbind(eye(4), zeros(4, 1))))
  initViQ <- matrix(inv(eye((5 * nQ)^2) - kron(BQ, BQ)) %*% as.vector(SQ), ncol = 5 * nQ, nrow = 5 * nQ)
  initViM <- diag(1 / diag(eye(size(BM, 1)) - BM ^ 2)) * SM
  
  # Output
  A <- blkdiag(A, BM, BQ)                # Measurement matrix
  Q <- blkdiag(Q, SM, SQ)                # Residual covariance matrix (transition)
  Z0 <- zeros(dim(A)[1], 1)              # States
  V0 <- blkdiag(V0, initViM, initViQ)    # Covariance of states
  
  
  return(list(A = A, C = C, Q = Q, R = R, Z0 = Z0, V0 = V0))
  
}