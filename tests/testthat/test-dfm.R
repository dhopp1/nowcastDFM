# test dataframe
test_data <- data.frame(
  date = as.Date(c("2006-12-01", "2007-01-01", "2007-02-01", "2007-03-01", "2007-04-01", "2007-05-01", "2007-06-01", "2007-07-01", "2007-08-01", "2007-09-01", "2007-10-01")),
  x = c(0.03861690, NA, NA, 0.02480291, NA, NA, 0.02795467, NA, NA, 0.04792785, NA),
  y = c(-0.041809423, 0.095962017, 0.177013620, -0.219734058, 0.139663588, 0.006332314, 0.035128206, 0.006955660, 0.025996546, -0.008288688, 0.020372225),
  z = c(0.012892695, -0.015530209, 0.013071097, 0.002588409, 0.051072682, -0.015054030, 0.001532301, 0.026254593, 0.013134543, 0.035182733, 0.018199641),
  aa = c(0.004666950, 0.009590594, -0.021090405, 0.044371479, -0.003273146, 0.018979616, 0.010749496, 0.004499497, 0.018575173, 0.005258734, 0.014724927)
)

old_data <- data.frame(test_data)
old_data[nrow(old_data), "y"] <- NA

output_dfm <- dfm(test_data)

test_that("DFM estimates", {
  expect_equal(output_dfm$convergence, 1)
})

test_that("predictions work", {
  expect_equal(
    round(predict_dfm(test_data, output_dfm, 3)[14,2], 3), 
    round(0.042574517, 3)
  )
})

test_that("gen_news works", {
  expect_equal(
    gen_news(old_data, test_data, output_dfm, "x", "2007-10-01")$impact_revisions, 
    0
  )
})


