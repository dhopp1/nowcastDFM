#' @title Viewing the impact of new data on a nowcast.
#' @description given and old and new dataset, will calculate the impact data releases and revisions have on the estimate of a target variable.
#' @param old_y dataframe of variables, size (n_obs, n_variables). Must include in 1st column a series of type date, called "date", all data already stationary.
#' @param new_y dataframe of variables, size (n_obs, n_variables). Must include in 1st column a series of type date, called "date", all data already stationary. Must contain same columns as old_y.
#' @param output_dfm list, the output of the \code{dfm()} function.
#' @param target_variable name of the target column.
#' @param target_period date of forecast to view impacts on.
#' @return A \code{list} containing the following elements:

#' \item{target_period}{same as input.}
#' \item{target_variable}{same as input.}
#' \item{y_old}{forecast for target variable with old data.}
#' \item{y_new}{forecast for target variable with new data.}
#' \item{forecast}{forecast of variables for target period. Only shows for variables that were newly published between old and new dataset.}
#' \item{actual}{actual published value of variables for target period. Only shows for variables that were newly published between old and new dataset.}
#' \item{weight}{weight of each data release}
#' \item{news_table}{table summarising forecast, actual, weight and impact of data releases}
#' \item{impact_revisions}{impact of data revisions on nowcast.}
#' \item{impact_releases}{impact of data releases on nowcast.}
#' \item{impact_total}{total impact (from data revision and data releases).}
#' 
#' @export

gen_news <- function(old_y, new_y, output_dfm, target_variable, target_period) {
  # making sure old data has same rows as new data
  old_y <- data.frame(date=new_y$date) %>% 
    left_join(old_y, by="date")
  
  data_old <- old_y
  data_new <- new_y
  
  ### Sort variables, first monthly variables, then quarterly variables
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
  for (i in 2:ncol(data_new)) {
    quarterly <- append(quarterly, is_quarterly(data_new[,1], data_new[,i]))
  }
  monthlies <- data_old[,which(quarterly == FALSE)]
  quarterlies <- data_old[,which(quarterly == TRUE)]
  column_names <- c(colnames(data_old)[which(quarterly == FALSE)], colnames(data_old)[which(quarterly == TRUE)])
  data_old <- cbind(monthlies, quarterlies)
  colnames(data_old) <- column_names
  
  monthlies <- data_new[,which(quarterly == FALSE)]
  quarterlies <- data_new[,which(quarterly == TRUE)]
  data_new <- cbind(monthlies, quarterlies)
  colnames(data_new) <- column_names
  
  t_nowcast <- which(data_new$date == target_period)
  
  # add 12 months to each dataset to allow for forecasting
  add_month <- function (X) {
    month <- as.numeric(substr(X, 6, 7))
    year <- as.numeric(substr(X, 1, 4))
    if (month == 12) {
      return (as.Date(paste0(year+1, "-01-01")))
    } else {
      return (as.Date(paste0(year, "-", month+1, "-01")))
    }
  }
  for (i in 1:12) {
    data_old[nrow(data_new) + 1, "date"] <- add_month(data_old[nrow(data_old), "date"])
    data_new[nrow(data_new) + 1, "date"] <- add_month(data_new[nrow(data_new), "date"])
  }
  
  # drop date column
  data_old <- data_old[,2:ncol(data_old)] 
  data_new <- data_new[,2:ncol(data_new)]
  i_series <- which(colnames(data_new) == target_variable)
  N <- ncol(data_new)
  
  # Update nowcast for target variable 'series' (i) at horizon 'target' (t)
  # Relate nowcast update into news from data releases:
  #   a. Compute the impact from data revisions
  #   b. Compute the impact from new data releases
  data_rev <- new_y
  data_rev[is.na(old_y)] <- NA
  
  # Compute news --------------------------------------------------------
  
  # Compute impact from data revisions
  results_old <- news_dfm(old_y, data_rev, output_dfm, target_variable, target_period)
  y_old <- results_old$y_old
  
  # Compute impact from data releases
  results_new <- news_dfm(data_rev, new_y, output_dfm, target_variable, target_period)
  y_rev <- results_new$y_old; y_new <- results_new$y_new
  actual <- results_new$actual; forecast <- results_new$fore; weight <- results_new$weight
  
  # Display output
  if (sum(is.na(forecast)) == length(forecast)) {
    print("No forecast was made")
    news_table <- NULL
    impact_revisions <- 0
    impact_releases <- 0
  } else {
    impact_revisions <- y_rev - y_old      # Impact from revisions
    news <- actual - forecast              # News from releases
    impact_releases <- sweep(weight, MARGIN = 1, news, "*")     # Impact of releases
    news_table <- data.frame(cbind(forecast, actual, weight, impact_releases), row.names = colnames(data_old))
    colnames(news_table) <- c("Forecast", "Actual", "Weight", "Impact")
    news_table[,"New Data"] <- as.numeric(as.logical(colSums(is.na(data_old) & !is.na(data_new))))
    impact_total <- impact_revisions + colSums(impact_releases, na.rm = T)
    
    print("Nowcast Impact Decomposition")
    print(paste("old nowcast: ", y_old * 100, "%", sep = ""))
    print(paste("new nowcast: ", y_new * 100, "%", sep = ""))
    print(paste("Impact from data revisions: ", sprintf("%.2f", impact_revisions * 100), "%", sep = ""))
    print(paste("Impact from data releases: ", 
                sprintf("%.2f", sum(news_table[, "Impact"] * 100, na.rm = TRUE)), "%", sep = ""))
    print(paste("Total impact: ", 
                sprintf("%.2f", (impact_revisions + sum(news_table[, "Impact"], na.rm = TRUE)) * 100),
                "%", sep = ""))
    print("Nowcast Detail Table")
    print(news_table[, c("Forecast", "Actual", "Weight", "Impact")])
  }
  
  return(list(target_period = target_period, target_variable = target_variable, y_old = y_old, y_new = y_new, forecast = forecast, actual = actual, weight = weight,
              news_table = news_table, impact_revisions = impact_revisions, 
              impact_releases = impact_releases,
              impact_total = impact_total))
  
}