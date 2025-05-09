test_that("rrda.fit returns expected structure (component = TRUE)", {
  set.seed(123)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  Y <- matrix(rnorm(100 * 10), nrow = 100)

  result <- rrda.fit(Y = Y, X = X, nrank = 1:3, lambda = c(0.1, 1), component = TRUE)

  expect_type(result, "list")
  expect_true("Bhat_comp" %in% names(result))
  expect_true("lambda" %in% names(result))
  expect_true("rank" %in% names(result))
  expect_true(is.list(result$Bhat_comp))
  expect_length(result$Bhat_comp, 2)  # 2 lambdas
})

test_that("rrda.fit returns matrix-form output (component = FALSE)", {
  set.seed(456)
  X <- matrix(rnorm(100 * 30), nrow = 100)
  Y <- matrix(rnorm(100 * 15), nrow = 100)

  result <- rrda.fit(Y = Y, X = X, nrank = c(1, 2), lambda = 1, component = FALSE)

  expect_type(result, "list")
  expect_true("Bhat_mat" %in% names(result))
  expect_length(result$Bhat_mat, 1)  # 1 lambda
  expect_true(is.list(result$Bhat_mat[[1]]))  # matrix list for each rank
})

test_that("rrda.fit stops with non-numeric nrank", {
  X <- matrix(rnorm(50), nrow = 10)
  Y <- matrix(rnorm(50), nrow = 10)

  expect_error(rrda.fit(Y, X, nrank = "a"))
})

test_that("rrda.fit stops with too large rank", {
  X <- matrix(rnorm(100), nrow = 10)
  Y <- matrix(rnorm(100), nrow = 10)

  expect_error(rrda.fit(Y, X, nrank = 20))
})

test_that("rrda.fit stops with negative lambda", {
  X <- matrix(rnorm(100), nrow = 10)
  Y <- matrix(rnorm(100), nrow = 10)

  expect_error(rrda.fit(Y, X, lambda = -1))
})
