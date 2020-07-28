news_dfm <- function(data_old, data_new, output_dfm, target_variable, target_period) {
  ###  This function calculates changes in news by using a given DFM results  structure.
  ###   It inputs two datasets, DFM parameters, target time index, and target variable index.
  ###   The function then produces Nowcast updates and decomposes the changes into news.
  ### 
  ###   Inputs:
  ###     data_old:  old data matrix (old vintage)
  ###     data_new:  new data matrix (new vintage)
  ###     output_dfm:  dfm.r output results
  ###     target_variable: name of target column
  ###     target_period: target date
  ###   
  ###  Outputs:
  ###       y_old:       old nowcast
  ###       y_new:       new nowcast
  ###       singlenews:  news for each data series
  ###       actual:      observed series release values
  ###       fore:        forecasted series values
  ###       weight:      news weight
  ###       t_miss:      time index for data releases
  ###       v_miss:      series index for data releases
  ###       innov:       difference between observed and predicted series values ("innovation")
  
  # making sure old data has same rows as new data
  data_old <- data.frame(date=data_new$date) %>% 
    left_join(data_old, by="date")
  t_nowcast <- which(data_new$date == target_period)
  # dropping date column
  data_old <- data_old[,2:ncol(data_old)]
  data_new <- data_new[,2:ncol(data_new)]
  i_series <- which(colnames(data_new) == target_variable)
  
  r <- dim(output_dfm$C)[2]
  N <- dim(data_new)[2]
  singlenews <- matrix(0, nrow = 1,  ncol = N)
  
  y_old <- y_new <- matrix(NA, ncol = length(i_series), nrow = 1)
  
  if (sum(is.na(data_new[t_nowcast, i_series])) == 0) {  ### No forecast case (values for i_series at time t_nowcast are already available)
    
    results_old <- kalman_filter_constparams(data_old, output_dfm, 0)
    for (i in 1:length(i_series)) {      # Loop for each target variable
      # (Observed value) - (predicted value)
      singlenews[, i_series[i]] <- data_new[t_nowcast, i_series[i]] - results_old$X_smooth[t_nowcast, i_series[i]]
      # Set predicted and observed y values
      y_old[1, i] <- results_old$X_smooth[t_nowcast, i_series[i]]
      y_new[1, i] <- data_new[t_nowcast, i_series[i]]
    }
    # Forecast-related output set to empty
    actual <- NULL; forecast <- NULL; weight <- NULL; t_miss <- NULL; v_miss <- NULL; innov <- NULL
    
  } else {   ### Forecast case (steps broken down into (A) and (B))
    
    # Initialize series mean/standard deviation respectively
    means <- output_dfm$means
    sdevs <- output_dfm$sdevs
    # Calculate indicators for missing values (1 if missing, 0 otherwise)
    miss_old <- is.na(data_old)
    miss_new <- is.na(data_new)
    # Indicator for missing--combine above information to single matrix where:
    # (i) -1: Value is in the old data, but missing in new data
    # (ii) 1: Value is in the new data, but missing in old data 
    # (iii) 0: Values are either both missing or both available in the datasets
    i_miss <- miss_old - miss_new
    t_miss <- which(i_miss == 1, arr.ind = TRUE)[, 1]   # Time/variable indices where case (ii) is true
    v_miss <- which(i_miss == 1, arr.ind = TRUE)[, 2]
    
    if (length(v_miss) == 0) {    ## Forecast subcase (A): No new information
      # Fill in missing variables using a Kalman filter
      results_old <- kalman_filter_constparams(data_old, output_dfm, 0)
      results_new <- kalman_filter_constparams(data_new, output_dfm, 0)
      # Set predicted and observed y values. New y value is set to old
      y_old <- results_old$X_smooth[t_nowcast, i_series]
      y_new <- y_old
      # y_new <- results_new$X_smooth[t_nowcast, i_series]
      # No news, so nothing returned for news-related output
      groupnews <- NULL; singlenews <- NULL; gain <- NULL; gainSer <- NULL
      actual <- NULL; forecast <- NULL; weight <- NULL; t_miss <- NULL; v_miss <- NULL; innov <- NULL
      
    } else {  ## Forecast subcase (b): new information
      # Difference between forecast time and new data time
      lag <- t_nowcast - t_miss
      # Gives biggest time interval between forecast and new data
      k <- max(abs(lag), (max(lag) - min(lag)))
      C <- output_dfm$C     # Measurement matrix
      R <- t(output_dfm$R)  # Covariance for measurement matrix residuals
      # Number of new events
      n_news <- length(lag)
      # Smooth old dataset
      results_old <- kalman_filter_constparams(data_old, output_dfm, k)
      Plag <- results_old$Plag
      # Smooth new dataset
      results_new <- kalman_filter_constparams(data_new, output_dfm, 0)
      # Subset for target variable and forecast time
      y_old <- results_old$X_smooth[t_nowcast, i_series]
      y_new <- results_new$X_smooth[t_nowcast, i_series]
      Vs <- results_old$Vsmooth[, , -1]
      P1 <- NULL  # Initialize projection onto updates
      
      # Cycle through total number of updates
      for (i in 1:n_news) {
        h <- abs(t_nowcast - t_miss[i])
        m <- max(t_miss[i], t_nowcast)
        # If location of update is later than the forecasting date
        if (t_miss[i] > t_nowcast) {
          Pp <- Plag[[h + 1]][, , m]  # P[1:r, h*r+1:h*r+r, m]'
        } else {
          Pp <- t(Plag[[h + 1]][, , m])  # P[1:r, h*r+1:h*r+r, m]
        }
        P1 <- cbind(P1, Pp %*% C[v_miss[i], 1:r])  # Projection on updates
      }
      
      innov <- matrix(NA, length(t_miss), 1)
      for (i in 1:length(t_miss)) {
        # Standardize predicted and observed values
        X_new_norm <- (data_new[t_miss[i], v_miss[i]] - means[v_miss[i]]) / sdevs[v_miss[i]]
        X_sm_norm <- (results_old$X_smooth[t_miss[i], v_miss[i]] - means[v_miss[i]]) / sdevs[v_miss[i]]
        # Innovation: gives [observed] data - [predicted data]
        innov[i] <- X_new_norm - X_sm_norm       
      }
      innov <- unlist(innov)
      
      ins <- dim(innov)[2]
      P2 <- NULL
      p2 <- NULL
      WW <- matrix(0, N, N)
      
      # Gives non-standardized series weights
      for (i in 1:length(lag)) {
        for (j in 1:length(lag)) {
          h <- abs(lag[i] - lag[j])
          m <- max(t_miss[i], t_miss[j])
          if (t_miss[j] > t_miss[i]) {
            Pp <- Plag[[h + 1]][, , m]  # P[1:r, h*r+1:(h+1)*r, m]'
          } else {
            Pp <- t(Plag[[h + 1]][, , m])  # P[1:r, h*r+1:(h+1)*r, m]
          }
          if (v_miss[i] == v_miss[j] && t_miss[i] != t_miss[j]) {
            WW[v_miss[i], v_miss[j]] <- 0;
          } else {
            WW[v_miss[i], v_miss[j]] <- R[v_miss[i], v_miss[j]]
          }
          p2 <- cbind(p2, C[v_miss[i], 1:r] %*% Pp %*% C[v_miss[j], 1:r] + WW[v_miss[i], v_miss[j]])
        }
        P2 <- rbind(P2, p2)
        p2 <- NULL
      }
      
      totnews <- matrix(NA, 1, length(i_series))
      temp <- array(NA, dim = c(1, n_news, length(i_series)))
      gain <- array(NA, dim = c(1, n_news, length(i_series)))
      for (i in 1:length(i_series)) {      # loop on v_news
        # Convert to real units (unstandardized data)
        totnews[1, i] <- unlist(sdevs[i_series[i]] * C[i_series[i], 1:r] %*% P1 %*% inv(P2) %*% innov)
        temp[1, , i] <- unlist(sdevs[i_series[i]]) %*% C[i_series[i], 1:r] %*% P1 %*% inv(P2) * innov
        gain[, , i] <- unlist(sdevs[i_series[i]]) * C[i_series[i], 1:r] %*% P1 %*% inv(P2)
      }
      
      # Initialize output objects
      singlenews <- array(NA, dim = c(max(t_miss) - min(t_miss) + 1, N, length(i_series)))
      actual     <- matrix(NA, N, 1)  # Actual forecasted values
      forecast   <- matrix(NA, N, 1)  # Forecasted values
      weight     <- array(NA, dim = c(N, 1, length(i_series)))
      # Fill in output values 
      for (i in 1:length(innov)) {
        actual[v_miss[i], 1] <- data_new[t_miss[i], v_miss[i]]  
        forecast[v_miss[i], 1] <- results_old$X_smooth[t_miss[i], v_miss[i]]
        for (j in 1:length(i_series)) {
          singlenews[t_miss[i] - min(t_miss) + 1, v_miss[i], j] <- temp[1, i, j]
          weight[v_miss[i], 1 , j] <- gain[, i, j] / unlist(sdevs[v_miss[i]])
        }
      }
      singlenews <- colSums(singlenews)  # Returns total news
      v_miss <- unique(v_miss)  
      
    }
    
  }
  
  return(list(y_old = y_old, y_new = y_new, singlenews = singlenews, actual = actual,
              fore = forecast, weight = weight, t_miss = t_miss, v_miss = v_miss, innov = innov))
  
}