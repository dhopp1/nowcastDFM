#' @title Estimating a dynamic factor model using the EM method.
#' @description Runs a DFM for the nowcast model on the transformed data at a certain data vintage. It uses an implementation through the EM algorithm. It relies on several functions to determine initial values, calculate the KF, the sequence of steps of the EM algorithm and criteria to determine convergence.
#' @param data matrix of variables, size (n_obs, n_variables). Must include in 1st column a series of type date, called "date", all data already stationary.
#' @param blocks Dataframe, size (n_variables, n_blocks). Note don't include date column in n_variables. Matrix of 1s or 0s for block loadings, i.e. 1 = included in block.
#' @param p number of lags in transition equation (AR element)
#' @param max_iter maximum number of iterations for EM (if no convergence)
#' @param threshold threshold for convergence of EM loop
#' @return A \code{list} containing the following elements:

#' \item{Xsmooth_std}{standardized Kalman-smoothed data where missing values are replaced by their expectation}
#' \item{Xsmooth}{Kalman-smoothed data where missing values are replaced by their expectation. In original input units.}
#' \item{Z}{smoothed states, rows give time, and columns are organized according to matrix C.}
#' \item{C}{measurement matrix, rows correspond to each series, and the columns are organized as, columns 1-20 give the factor loadings. For example, 1-5 give loadings for the first, and are organized in reverse-chronological order (f^G_t, f^G_t-1, f^G_t-2, f^G_t-3, f^G_t-4), Columns 6-10, 11-15, and 16-20 give loadings for the second, third, and fourth blocks respectively.}
#' \item{R}{covariance for measurement matrix residuals.}
#' \item{A}{transition matrix, a square matrix that follows the same organization scheme as matrix C's columns. Identity matrices are used to account for matching terms on the left and righthand side. For example, we place an I4 matrix to account for matching (f_t-1; f_t-2; f_t-3; f_t-4) terms.}
#' \item{Q}{covariance for transition equation residuals}
#' \item{means}{means of each column.}
#' \item{sdevs}{standard deviations of each column.}
#' \item{Z0}{initial value of state.}
#' \item{V0}{initial value of covariance matrix}
#' \item{p}{number of lags in transition equation (AR element).}
#' \item{model}{names of features input to the model.}
#' \item{blocks}{same as parameter passed in.}
#' \item{num_vars}{number of features estimated in the model.}
#' \item{num_iter}{number of iterations for log likelihood to converge or hit maximum.}
#' \item{convergence}{1 if algorithm converged successfully (given max_iter).}
#' \item{loglik}{log likelihood of last iteration.}
#' \item{LL}{sequence of log likelihoods per iteration.}
#' \item{data}{data passed to the model.}
#' 
#' @import matlab
#' 
#' @export

dfm <- function(data, blocks, p, max_iter=5000, threshold=1e-5) {
  ### Sort variables, first monthly variables, then quarterly variables
  orig_data <- data
  
  is_quarterly <- function(dates, series) {
    tmp <- data.frame(dates, series) %>% 
      dplyr::filter(!is.na(series)) %>% 
      select(dates) %>% pull
    if (identical((sapply(tmp, function(x) substr(x, 6, 7)) %>% unique %>% sort), c("03", "06", "09", "12"))) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  }
  quarterly <- c(FALSE)
  for (i in 2:ncol(data)) {
    quarterly <- append(quarterly, is_quarterly(data[,1], data[,i]))
  }
  monthlies <- data[,which(quarterly == FALSE)]
  quarterlies <- data[,which(quarterly == TRUE)]
  data <- cbind(monthlies, quarterlies)
  # drop date column
  data <- data[,2:ncol(data)]
  quarterly <- quarterly[2:length(quarterly)]
  index_freq <- as.integer(!quarterly)
  
  
  ### Obtain the characteristics of the model from the catalogue and the data
  num_obs <- nrow(data)
  nM <- sum(!quarterly)
  nQ <- sum(quarterly)
  num_blocks <- ncol(blocks)
  R_mat = matrix(c(2, 3, 2, 1, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1), ncol = 5)   # Quarterly-monthly aggregation scheme
  q <- zeros(4, 1)
  
  ### Standardize data
  means <- data %>% summarise_all(mean, na.rm = T)
  sdevs <- data %>% summarise_all(sd, na.rm = T)
  data_std <- data %>% mutate_all(~ scale(.))
  index_na <- is.na(data)
  
  ### Calculate initial values
  init <- init_conds(data_std, p, blocks, R_mat, q, nM, nQ, index_freq)
  A <- init$A; C <- init$C; Q <- init$Q; R <- init$R; Z0 <- init$Z0; V0 <- init$V0
  
  ## Initialize EM loop values
  prev_loglik <- -1e6
  num_iter <- 0
  LL <- -1e6
  converged <- 0
  y <- t(data_std)   ## y for the estimation is WITH missing data
  
  ## EM LOOP ----------------------------------------------------------------
  
  # The model can be written as
  # y = C*Z + e;
  # Z = A*Z(-1) + v
  # where y is NxT, Z is (pr)xT, etc.
  
  # Remove the leading and ending NAs
  y_est <- data_std %>%
    mutate(na_row = ifelse(rowSums(is.na(.)) == ncol(.), 1, 0)) %>%
    mutate(row_num = 1:nrow(.), row_num_inv = nrow(.):1) %>%
    mutate(beginning = ifelse(cumsum(na_row) == row_num, 1, 0)) %>%
    mutate(ending = ifelse(rev(cumsum(rev(na_row))) == row_num_inv, 1, 0)) %>%
    dplyr::filter(beginning != 1 & ending != 1) %>%
    select(-(na_row:ending)) %>%
    t(.)
  
  while(!converged & num_iter <= max_iter) {
    
    em_output <- EM_step(y_est, A, C, Q, R, Z0, V0, p, blocks, R_mat, q, nM, nQ, index_freq)
    A <- em_output$A_new; C <- em_output$C_new; Q <- em_output$Q_new; R <- em_output$R_new
    Z0 <- em_output$Z0; V0 <- em_output$V0; loglik <- em_output$loglik
    
    em_conv <- EM_convergence(loglik, prev_loglik, threshold)
    converged <- em_conv$converged
    
    if ((mod(num_iter, 20) == 0) & (num_iter > 0)) {   # Print a message every 20 iteratios
      message("Now running iteration number ", num_iter, " out of a maximum of ", max_iter)
      message("Loglik: ", sprintf("%.4f", loglik), "; % change: ", 
              sprintf("%.4f", 100 * (loglik - prev_loglik)/prev_loglik), "%")
    }
    
    LL <- c(LL, loglik)
    prev_loglik <- loglik
    num_iter <- num_iter + 1
    
  }
  
  # Final run of the Kalman filter
  kf_output <- kalman_filter(y, A, C, Q, R, Z0, V0)
  Zsmooth <- t(kf_output$Zsmooth)
  Xsmooth_std <- Zsmooth[2:nrow(Zsmooth), ] %*% t(C)
  
  # Create list with the results
  nowcast <- list()
  nowcast$Xsmooth_std <- Xsmooth_std
  nowcast$Xsmooth <- repmat(as.numeric(sdevs), num_obs, 1) * Xsmooth_std + repmat(as.numeric(means), num_obs, 1)
  nowcast$Z <- Zsmooth[2:nrow(Zsmooth), ]
  nowcast$C <- C
  nowcast$R <- R
  nowcast$A <- A
  nowcast$Q <- Q
  nowcast$means <- means
  nowcast$sdevs <- sdevs
  nowcast$Z0 <- Z0
  nowcast$V0 <- V0
  nowcast$p <- p
  nowcast$model <- colnames(data)[1:length(colnames(data))]
  nowcast$blocks <- blocks
  nowcast$num_vars <- nM + nQ
  nowcast$num_iter <- num_iter
  nowcast$convergence <- converged
  nowcast$loglik <- loglik
  nowcast$LL <- LL[2:length(LL)]
  nowcast$data <- orig_data
  
  return(nowcast)
  
}


#' @title Predictions from an estimated dynamic factor model.
#' @description runs a dataset through a previously estimated DFM to obtain predictions for all missing values in the series.
#' @param data matrix of variables, size (n_obs, n_variables). Must include in 1st column a series of type date, called "date", all data already stationary.
#' @param output_dfm list, the output of the \code{dfm()} function.
#' @param months_ahead number of months ahead to forecast.
#' @param lag number of lags for the kalman filter
#' @return dataframe with all missing values filled + predictions.
#' 
#' @export

predict_dfm <- function(data, output_dfm, months_ahead=3, lag=0) {
  # add months
  add_month <- function (X) {
    month <- as.numeric(substr(X, 6, 7))
    year <- as.numeric(substr(X, 1, 4))
    if (month == 12) {
      return (as.Date(paste0(year+1, "-01-01")))
    } else {
      return (as.Date(paste0(year, "-", month+1, "-01")))
    }
  }
  output <- data
  for (i in 1:months_ahead) {
    output[nrow(output) + 1, "date"] <- add_month(output[nrow(output), "date"])
  }
  dates <- output[,"date"]
  output <- output[,2:ncol(output)]
  
  # prediction with constant parameters
  preds <- kalman_filter_constparams(output, output_dfm, lag=lag)$X_smooth %>% 
    data.frame
  preds$date <- dates
  preds <- preds[,c(ncol(preds), 1:(ncol(preds)-1))]
  colnames(preds) <- colnames(data)
  
  return (preds)
}