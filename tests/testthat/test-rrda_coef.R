test_that("rrda.coef returns Bhat matrix for valid input", {
  set.seed(123)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  Y <- matrix(rnorm(100 * 10), nrow = 100)

  # Fit RRDA model
  fit_result <- rrda.fit(Y = Y, X = X, nrank = 1:3, lambda = c(0.1, 1), component = TRUE)

  # Convert to matrix
  Bhat_mat <- rrda.coef(Bhat = fit_result, nrank = 3, lambda = 0.1)

  expect_type(Bhat_mat, "list")
  expect_true("rank3" %in% names(Bhat_mat[[1]]))
  expect_true(is.matrix(Bhat_mat[[1]]$rank3))
})

test_that("rrda.coef throws error for invalid nrank", {
  set.seed(123)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  Y <- matrix(rnorm(100 * 10), nrow = 100)

  fit_result <- rrda.fit(Y = Y, X = X, nrank = 1:3, lambda = 0.1, component = TRUE)

  expect_error(rrda.coef(Bhat = fit_result, nrank = 999))  # 無効なrank
})

test_that("rrda.coef gives warning for partially valid nrank", {
  set.seed(456)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  Y <- matrix(rnorm(100 * 10), nrow = 100)

  fit_result <- rrda.fit(Y = Y, X = X, nrank = 1:3, lambda = 0.1, component = TRUE)

  expect_warning({
    Bhat_mat <- rrda.coef(Bhat = fit_result, nrank = c(2, 999))
  })
  expect_true("rank2" %in% names(Bhat_mat[[1]]))
})

test_that("rrda.coef handles input with Bhat_mat", {
  set.seed(789)
  X <- matrix(rnorm(100 * 15), nrow = 100)
  Y <- matrix(rnorm(100 * 5), nrow = 100)

  fit_result <- rrda.fit(Y = Y, X = X, nrank = 1:2, lambda = 0.5, component = FALSE)  # matrix形式で保存

  coef_result <- rrda.coef(Bhat = fit_result, nrank = 2, lambda = 0.5)

  expect_type(coef_result, "list")
  expect_true("rank2" %in% names(coef_result[[1]]))
})

