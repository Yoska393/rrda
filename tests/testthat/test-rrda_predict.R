test_that("rrda.predict returns predicted matrix in expected structure", {
  set.seed(123)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  Y <- matrix(rnorm(100 * 10), nrow = 100)

  # Fit the model
  Bhat <- rrda.fit(Y = Y, X = X, nrank = 1:3, lambda = c(0.1, 1))

  # Predict
  Yhat <- rrda.predict(Bhat = Bhat, X = X, nrank = 3, lambda = 1)

  expect_type(Yhat, "list")
  expect_true("Yhat_mat" %in% names(Yhat))
  expect_true(is.list(Yhat$Yhat_mat))
  expect_true("rank3" %in% names(Yhat$Yhat_mat[[1]]))
  expect_true(is.matrix(Yhat$Yhat_mat[[1]]$rank3))
  expect_equal(dim(Yhat$Yhat_mat[[1]]$rank3), dim(Y))
})

test_that("rrda.predict throws error if Bhat_mat is passed", {
  set.seed(456)
  X <- matrix(rnorm(100 * 15), nrow = 100)
  Y <- matrix(rnorm(100 * 5), nrow = 100)

  # component = FALSE → Bhat_mat に保存される
  Bhat <- rrda.fit(Y = Y, X = X, nrank = 1, lambda = 0.5, component = FALSE)

  expect_error(rrda.predict(Bhat = Bhat, X = X, nrank = 1, lambda = 0.5),
               "Bhat_mat is detected")
})

test_that("rrda.predict gives warning for partially valid nrank or lambda", {
  set.seed(789)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  Y <- matrix(rnorm(100 * 5), nrow = 100)

  Bhat <- rrda.fit(Y = Y, X = X, nrank = 1:2, lambda = c(0.1, 1))

  expect_warning({
    Yhat <- rrda.predict(Bhat = Bhat, X = X, nrank = c(2), lambda = c(0.1, 9))
  })
  expect_warning({
    Yhat <- rrda.predict(Bhat = Bhat, X = X, nrank = c(2, 999), lambda = c(0.1))
  })
  expect_true("rank2" %in% names(Yhat$Yhat_mat[[1]]))
})

test_that("rrda.predict throws error for invalid nrank or lambda", {
  set.seed(999)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  Y <- matrix(rnorm(100 * 10), nrow = 100)

  Bhat <- rrda.fit(Y = Y, X = X, nrank = 1:2, lambda = c(0.1, 1))

  expect_error(rrda.predict(Bhat = Bhat, X = X, nrank = 999), "nrank must be included")
  expect_error(rrda.predict(Bhat = Bhat, X = X, lambda = 999), "lambda must be included")
})

test_that("rrda.predict returns correctly scaled-back Yhat", {
  set.seed(321)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  Y <- matrix(rnorm(100 * 5), nrow = 100)

  Bhat <- rrda.fit(Y = Y, X = X, nrank = 2, lambda = 0.1,
                   center.X = TRUE, center.Y = TRUE,
                   scale.X = TRUE, scale.Y = TRUE)

  Yhat <- rrda.predict(Bhat = Bhat, X = X, nrank = 2, lambda = 0.1)
  yhat_mat <- Yhat$Yhat_mat[[1]]$rank2

  expect_equal(dim(yhat_mat), dim(Y))
  expect_true(is.numeric(yhat_mat))
})
