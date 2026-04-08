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

test_that("rrda.top returns expected structure and classes", {
  set.seed(123)
  sim <- rdasim1(n = 100, p = 50, q = 20, k = 5)
  X <- sim$X
  Y <- sim$Y

  res <- rrda.top(Y = Y, X = X, nrank = 3, lambda = 1, mx = 10, my = 8, title = "custom")
  expect_type(res, "list")
  expect_true(all(c("heatmap","B_sub","top_x","top_y","b1_sub","b2_sub") %in% names(res)))

  # heatmap object
  expect_s3_class(res$heatmap, "pheatmap")

  # Submatrix dimensions
  expect_equal(dim(res$B_sub), c(10, 8))

  # top_x / top_y lengths
  expect_length(res$top_x, 10)
  expect_length(res$top_y, 8)

  # Consistency of row names
  expect_setequal(rownames(res$b1_sub), res$top_x)
  expect_setequal(rownames(res$b2_sub), res$top_y)
})

test_that("rrda.top clamps mx/my < 1 to 1 without error", {
  set.seed(123)
  sim <- rdasim1(n = 50, p = 12, q = 7, k = 3)
  X <- sim$X; Y <- sim$Y

  res <- rrda.top(Y, X, nrank = 2, mx = 0, my = 0)
  expect_equal(dim(res$B_sub), c(1, 1))
  expect_length(res$top_x, 1)
  expect_length(res$top_y, 1)
})

test_that("rrda.top warns and caps mx/my when exceeding ncol(X)/ncol(Y)", {
  set.seed(123)
  sim <- rdasim1(n = 60, p = 15, q = 9, k = 3)
  X <- sim$X; Y <- sim$Y

  expect_warning(
    res <- rrda.top(Y, X, nrank = 2, mx = ncol(X) + 5, my = ncol(Y) + 10),
    "mx is larger|my is larger"
  )
  expect_equal(dim(res$B_sub), c(ncol(X), ncol(Y)))
  expect_length(res$top_x, ncol(X))
  expect_length(res$top_y, ncol(Y))
})

test_that("rrda.top sets default nrank and emits message", {
  set.seed(123)
  sim <- rdasim1(n = 40, p = 10, q = 6, k = 2)
  X <- sim$X; Y <- sim$Y

  expect_message(
    res <- rrda.top(Y, X, nrank = NULL, lambda = 1, mx = 3, my = 3),
    "nrank is set to default"
  )
  # Default should be min(5, min(dim(X), dim(Y)))
  expect_equal(dim(res$B_sub), c(3, 3))
})

test_that("rrda.top validates arguments and throws informative errors", {
  set.seed(123)
  sim <- rdasim1(n = 30, p = 8, q = 5, k = 2)
  X <- sim$X; Y <- sim$Y

  # Non-numeric nrank
  expect_error(rrda.top(Y, X, nrank = "a", lambda = 1), "nrank must be numeric")

  # rank too large
  too_big <- min(dim(X), dim(Y)) + 1
  expect_error(rrda.top(Y, X, nrank = too_big, lambda = 1), "rank\\(B\\) must be less than or equal")

  # Negative lambda
  expect_error(rrda.top(Y, X, nrank = 1, lambda = -0.1), "non-negative")
})

test_that("rrda.top is deterministic with fixed seed for top indices size", {
  set.seed(2024)
  sim <- rdasim1(n = 120, p = 40, q = 18, k = 4)
  X <- sim$X; Y <- sim$Y

  res1 <- rrda.top(Y, X, nrank = 2, lambda = 1, mx = 7, my = 9)
  set.seed(2024)
  res2 <- rrda.top(Y, X, nrank = 2, lambda = 1, mx = 7, my = 9)

  # Check consistency of lengths and dimensions (not exact values)
  expect_equal(length(res1$top_x), length(res2$top_x))
  expect_equal(length(res1$top_y), length(res2$top_y))
  expect_equal(dim(res1$B_sub), dim(res2$B_sub))
})




