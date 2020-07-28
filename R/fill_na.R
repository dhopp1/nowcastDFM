#' @importFrom signal filter
#' @importFrom pracma cubicspline

fill_na <- function(X) {
  ###
  ### Function that fills missing values from the data. Useful for calculating starting values through
  ### principal components, who can't take missing values.
  ###
  ### Inputs:
  ###    - X: standardised data (required, no default)
  ###
  ### Outputs:
  ###    - data_full: dataset with missing values replaced with cubic spline (between observed values) and 
  ###                 1-D digital filter fitted values for missing tails calculated at the variable level.
  ###                 Rows at tails with more than 80% missing were dropped.
  ###    - ind_na: indicator matrix of NAs
  
  k <- 3
  
  temp <- X %>%
    mutate(na_row = ifelse(rowSums(is.na(.)) > 0.8 * ncol(X), 1, 0)) %>%
    mutate(row_num = 1:nrow(.), row_num_inv = nrow(.):1) %>%
    mutate(beginning = ifelse(cumsum(na_row) == row_num, 1, 0)) %>%
    mutate(ending = ifelse(rev(cumsum(rev(na_row))) == row_num_inv, 1, 0)) %>%
    dplyr::filter(beginning != 1 & ending != 1) %>%
    select(-(na_row:ending))
  ind_na <- is.na(temp)
  
  for (i in 1:ncol(temp)) {
    tempi <- temp[, i]
    ind_na_i <- is.na(tempi)
    t1 <- min(which(!is.na(tempi)))
    t2 <- max(which(!is.na(tempi)))
    tempi[t1:t2] <- cubicspline(x = which(!is.na(tempi)), y = tempi[which(!is.na(tempi))], xi = t1:t2)
    ind_na_i <- is.na(tempi)
    tempi[ind_na_i] = median(tempi, na.rm = T)
    tempi_MA <- signal::filter(filt = ones(2 * k + 1, 1) / (2 * k + 1), a = 1, 
                       x = c(tempi[1] * ones(k, 1), tempi, tempi[nrow(tempi)] * ones(k, 1)))
    tempi_MA <- tempi_MA[(2 * k + 1):length(tempi_MA)]
    tempi[ind_na_i] <- tempi_MA[ind_na_i] 
    temp[, i] <- tempi
  }
  
  return(list(data_full = temp, ind_na = ind_na))
  
}