# nowcastDFM
Dynamic factor models for R. Adapted from [Bok et al. 2017](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr830.pdf), [code](https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting).

# Installation
`devtools::install_github("dhopp1-UNCTAD/nowcastDFM")`

# Functionality
- `dfm`: estimate a dynamic factor model using the EM method. `?dfm` for more info.
- `predict_dfm`: obtain predictions from a previously estimated model. `?predict_dfm` for more info.
- `gen_news`: obtain impacts of new data releases and revisions on the forecast of a target variable. `?gen_news` for more info.
