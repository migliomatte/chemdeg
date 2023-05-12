# FOMT
test_that("FOMT works with dtataframe with 3 columns(error)", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  expect_no_error(FOMT <- FOMT(dframe))
  regr <- nls(conc ~ 1 - (1 - exp(-k * t))^n, start = par_est_FOMT(t, conc),
              weights = 1 / err^2)
  expect_equal(FOMT$m$getAllPars(), regr$m$getAllPars(), ignore_attr = TRUE)
})

test_that("FOMT works with dtataframe with 2 columns (no error)", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  dframe <- data.frame(t, conc)
  expect_no_error(FOMT <- FOMT(dframe))
  regr <- nls(conc ~ 1 - (1 - exp(-k * t))^n, start = par_est_FOMT(t, conc))
  expect_equal(FOMT$m$getAllPars(), regr$m$getAllPars(), ignore_attr = TRUE)
})

test_that("FOMT works with dtataframe with 2 columns (no error)", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err, err)
  expect_error(FOMT(t))
  expect_error(FOMT(dframe))
})

# FOMTm
test_that("FOMTm resturns correct values", {
  k <- 0.5
  t <- c(0, 1, 2)
  n <- 3
  expect_equal(FOMTm(k, t, n), c(1, 0.939, 0.747), tolerance = 0.001)
})

test_that("FOMTm can be used as formula in nls function", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  expect_no_error(fit <- nls(conc ~ FOMTm(k, t, n),
                             start = list(k = 0.5, n = 599)))
})


# par_est_FOMT
test_that("par_est_FOMT works with a data frame with 2 columns", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  dframe <- data.frame(t, conc)
  expect_no_error(par_est_FOMT(dframe))
  expect_equal(par_est_FOMT(dframe), c(k = 0.34657, n = 30.72000),
               tolerance = 1e-4)
})

test_that("par_est_FOMT works with a two arrays", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  expect_no_error(par_est_FOMT(t, conc))
  expect_equal(par_est_FOMT(t, conc), c(k = 0.34657, n = 30.72000),
               tolerance = 1e-4)
})

test_that("par_est_FOMT returns the same results with dataframe and arrays", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
  dframe <- data.frame(t, conc)
  expect_no_error(par_est_FOMT(dframe))
  expect_no_error(par_est_FOMT(t, conc))
  expect_equal(par_est_FOMT(dframe), par_est_FOMT(t, conc))
})
