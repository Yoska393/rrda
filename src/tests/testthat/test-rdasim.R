test_that("rdasim1 generates data with correct dimensions", {
  set.seed(123)
  sim <- rdasim1(n = 10, p = 5, q = 3, k = 2)

  expect_type(sim, "list")
  expect_true(all(c("X", "Y", "H", "theta.y") %in% names(sim)))
  expect_equal(dim(sim$X), c(10, 5))
  expect_equal(dim(sim$Y), c(10, 3))
  expect_equal(dim(sim$H), c(10, 2))
  expect_equal(dim(sim$theta.y), c(2, 3))
})

test_that("rdasim2 generates data with correct dimensions", {
  set.seed(123)
  sim <- rdasim2(n = 20, p = 8, q = 4, k = 3)

  expect_type(sim, "list")
  expect_true(all(c("X", "Y", "B", "E") %in% names(sim)))
  expect_equal(dim(sim$X), c(20, 8))
  expect_equal(dim(sim$Y), c(20, 4))
  expect_equal(dim(sim$B), c(8, 4))
  expect_equal(dim(sim$E), c(20, 4))
})


test_that("rdasim2 works with non-zero off-diagonal covariance", {
  set.seed(456)
  sim <- rdasim2(n = 30, p = 6, q = 3, k = 2, xofd = 0.5)

  expect_type(sim$X, "double")
  expect_equal(dim(sim$X), c(30, 6))
})
