## code to prepare `fomtdata` dataset goes here

fomtdata <- data.frame(
  time_h = c(0, 5, 24, 48, 72, 96, 144, 168),
  tCQA_AA = c(1, 0.997, 0.994, 0.424, 0.190, 0.072, 0.053, 0.001)
)

usethis::use_data(fomtdata, overwrite = TRUE)
