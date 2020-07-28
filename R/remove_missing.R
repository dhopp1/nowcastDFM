remove_missing <- function(yt, C, R) {
  
  #   This function eliminates the rows in y and matrices C and R that correspond to missing data in y
  #
  # Input:
  #   - yt: vector of observations at time t
  #   - C:  measurement matrix
  #   - R:  covariance for measurement matrix residuals
  #
  # Output:
  #   - yt: vector of observations at time t (reduced)     
  #   - C:  measurement matrix (reduced)     
  #   - R:  covariance for measurement matrix residuals
  #   - L:  used to restore standard dimensions(n x w) where w is the number of available data in y
  
  # Returns 1 for nonmissing series
  ix <- !is.na(yt)
  
  # Index for columns with nonmissing variables
  e <- eye(length(yt))
  L <- e[, ix]
  
  # Removes missing series
  yt <- yt[ix]

  # Removes missing series from observation matrix
  C <- C[ix, ]
  
  # Removes missing series from transition matrix
  R <- R[ix, ix]
  
  # Prepare output
  return(list(yt = yt, C = C, R = R, L = L))
  
}