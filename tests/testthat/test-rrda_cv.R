test_that("rrda.cv returns expected structure", {
  set.seed(123)
  X <- matrix(rnorm(100 * 50), nrow = 100)
  Y <- matrix(rnorm(100 * 20), nrow = 100)

  result <- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 3, verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(c("MSE", "SEM", "rank", "lambda", "opt_min", "opt_lambda.1se", "opt_rank.1se") %in% names(result)))
  expect_true(is.data.frame(result$MSE))
  expect_true(is.data.frame(result$SEM))
  expect_true(is.numeric(result$rank))
  expect_true(is.numeric(result$lambda))
})

test_that("rrda.cv returns proper optimal values", {
  set.seed(456)
  X <- matrix(rnorm(60 * 30), nrow = 60)
  Y <- matrix(rnorm(60 * 10), nrow = 60)

  result <- rrda.cv(Y = Y, X = X, maxrank = 3, nfold = 2, verbose = FALSE)

  expect_true(is.list(result$opt_min))
  expect_true(is.numeric(result$opt_min$rank))
  expect_true(is.numeric(result$opt_min$lambda))
})

test_that("rrda.cv throws error with negative lambda", {
  set.seed(789)
  X <- matrix(rnorm(50 * 20), nrow = 50)
  Y <- matrix(rnorm(50 * 10), nrow = 50)

  expect_error(rrda.cv(Y = Y, X = X, lambda = c(-0.1, 0.2), maxrank = 2, nfold = 2, verbose = FALSE))
})

test_that("rrda.cv handles rank too large with warning", {
  X <- matrix(rnorm(50 * 10), nrow = 50)
  Y <- matrix(rnorm(50 * 8), nrow = 50)

  expect_warning({
    result <- rrda.cv(Y = Y, X = X, maxrank = 100, nfold = 2, verbose = FALSE)
  })
  expect_true(result$opt_min$rank <= min(dim(X), dim(Y)))
})

test_that("rrda.cv throws error for nfold < 2", {
  sim <- rdasim1(n = 30, p = 10, q = 5, k = 2)
  X <- sim$X
  Y <- sim$Y

  expect_error(rrda.cv(Y = Y, X = X, nfold = 1),
               "nfold must be at least 2")
})

