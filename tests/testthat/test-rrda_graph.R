test_that("rrda.plot returns a ggplot object with default parameters", {
  set.seed(123)
  sim <- rdasim1(n = 100, p = 50, q = 20, k = 5)
  X <- sim$X
  Y <- sim$Y

  cv_result <- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 3, verbose = FALSE)

  p <- rrda.plot(cv_result)

  expect_s3_class(p, "ggplot")
})

test_that("rrda.plot works with custom rank and lambda range", {
  set.seed(123)
  sim <- rdasim1(n = 100, p = 50, q = 20, k = 5)
  X <- sim$X
  Y <- sim$Y

  cv_result <- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 2, verbose = FALSE)

  lambda_range <- range(cv_result$lambda)
  p <- rrda.plot(cv_result, nrank = c(2, 4), min_l = lambda_range[1], max_l = lambda_range[2], show_error_bar = TRUE)

  expect_s3_class(p, "ggplot")
})

test_that("rrda.heatmap returns a ggplot object with default parameters", {
  set.seed(123)
  sim <- rdasim1(n = 100, p = 50, q = 20, k = 5)
  X <- sim$X
  Y <- sim$Y

  cv_result <- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 2, verbose = FALSE)

  p <- rrda.heatmap(cv_result)

  expect_s3_class(p, "ggplot")
})

test_that("rrda.heatmap returns heatmap with highlighted points", {
  set.seed(123)
  sim <- rdasim1(n = 100, p = 50, q = 20, k = 5)
  X <- sim$X
  Y <- sim$Y

  cv_result <- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 2, verbose = FALSE)

  p <- rrda.heatmap(cv_result, highlight_min = TRUE)

  expect_s3_class(p, "ggplot")
})
