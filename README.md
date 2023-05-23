
# nowcastDFM
Run dynamic factor models (DFM) in R. Adapted from [Bok et al. 2017](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr830.pdf), [MATLAB code](https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting). The package provides the ability to estimate a DFM model using the expectation–maximization method, obtain predictions from estimated models, and obtain the impact of new data releases on model predictions. On [CRAN](https://cloud.r-project.org/web/packages/nowcastDFM/index.html).

# Installation
```R
install.packages("nowcastDFM")`
```
If this does not work, you can install directly from Github with: 
```R
install.packages("devtools")
devtools::install_github("dhopp1/nowcastDFM")
```

# Functionality
- `dfm`: estimate a dynamic factor model using the EM method. `?dfm` for more info.
- `predict_dfm`: obtain predictions from a previously estimated model. `?predict_dfm` for more info.
- `gen_news`: obtain impacts of new data releases and revisions on the forecast of a target variable. `?gen_news` for more info.

# Example
Given `data` is a dataframe (not a tibble) with a `date` column  and 4 columns for various seasonally adjusted growth rates of economic series with missing values of `NA`:
```R
library(nowcastDFM)

# estimate a DFM with one block for all variables
output_dfm <- dfm(data) 

# estimate a DFM with two different blocks
blocks <- data.frame(block_1 = c(1,1,1,0), block_2 = c(0,0,1,1)) # defining two blocks
output_dfm <- dfm(data, blocks = blocks)

# get predictions from estimated DFM for the following 3 months
# new data is dataframe with same columns as data the model was trained on, but newer data
predictions <- predict_dfm(new_data, output_dfm, months_ahead = 3)

# get impact of new data on predictions for a particular variable and time period
# old_data and new_data are dataframes with same columns as the data the model was trained on, but with older and newer data
news <- gen_news(old_data, new_data, output_dfm, target_variable = "target_name", target_period = "2020-01-01")
```
