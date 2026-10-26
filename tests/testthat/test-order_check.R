sim_data <- function(R, n, seed) {
  set.seed(seed)
  X <- as.data.frame(matrix(stats::rnorm(n * nrow(R)), n) %*% chol(unclass(R)))
  names(X) <- rownames(R)
  X
}

test_that("order_check finds both minima and the order effect", {
  skip_on_cran()
  R <- two_trait_cor(phi = 0.1, lambda_A = 0.68, lambda_B = 0.72)
  X <- sim_data(R, 2000, 1)
  oc <- order_check(X, blocks = list(A = paste0("A", 1:4), B = paste0("B", 1:4)), nstart = 10)
  expect_s3_class(oc, "order_check")
  expect_equal(nrow(oc$solutions), 2)
  expect_equal(oc$solutions$follows[1], "B")
  expect_true(oc$order_dependent)
  expect_equal(oc$default_first_A, 2L)
  expect_output(print(oc), "Item order changes the solution")
})

test_that("order_check reports one solution when the traits are highly correlated", {
  skip_on_cran()
  R <- two_trait_cor(phi = 0.8, lambda_A = 0.7, m_A = 4)
  X <- sim_data(R, 1000, 2)
  oc <- order_check(X, blocks = list(dep = paste0("A", 1:4), anx = paste0("B", 1:4)), nstart = 5)
  expect_equal(nrow(oc$solutions), 1)
  expect_false(oc$order_dependent)
  expect_output(print(oc), "share_dep")
})

test_that("order_check runs with WLSMV on ordinal items", {
  skip_on_cran()
  R <- two_trait_cor(phi = 0.3, lambda_A = 0.7, m_A = 4)
  X <- sim_data(R, 800, 3)
  X[] <- lapply(X, function(v) cut(v, c(-Inf, -0.8, 0, 0.8, Inf), labels = FALSE))
  oc <- order_check(X, blocks = list(paste0("A", 1:4), paste0("B", 1:4)),
                    estimator = "WLSMV", ordered = TRUE, nstart = 3)
  expect_gte(nrow(oc$solutions), 1)
  expect_true("cfi_robust" %in% names(oc$solutions))
})

test_that("order_check checks its input", {
  X <- data.frame(a = 1:5, b = 1:5)
  expect_error(order_check(X, list("a")), "two character vectors")
  expect_error(order_check(X, list("a", "z")), "not found")
})

test_that("ULS fits correlations, so rescaling an item changes nothing", {
  skip_on_cran()
  X <- sim_data(two_trait_cor(phi = 0.2, lambda_A = 0.7, m_A = 4), 2000, 4)
  bl <- list(A = paste0("A", 1:4), B = paste0("B", 1:4))
  a <- order_check(X, bl, estimator = "ULS", nstart = 3)
  X$A1 <- 3 * X$A1
  b <- order_check(X, bl, estimator = "ULS", nstart = 3)
  expect_equal(a$solutions$share_A, b$solutions$share_A, tolerance = 1e-6)
  m <- order_check(cov(X), bl, estimator = "ULS", nstart = 3, nobs = nrow(X))
  expect_equal(m$solutions$share_A, b$solutions$share_A, tolerance = 1e-6)
})

test_that("a covariance matrix gives the same solutions as raw data", {
  skip_on_cran()
  X <- sim_data(two_trait_cor(phi = 0.1, lambda_A = 0.68, lambda_B = 0.72), 2000, 1)
  bl <- list(A = paste0("A", 1:4), B = paste0("B", 1:4))
  d <- order_check(X, bl, nstart = 5)
  m <- order_check(cov(X), bl, nstart = 5, nobs = nrow(X))
  expect_equal(m$solutions$fmin, d$solutions$fmin, tolerance = 1e-8)
  expect_equal(m$default_first_A, d$default_first_A)
  expect_true(all(d$solutions$min_uniqueness > 0))
  expect_error(order_check(cov(X), bl), "nobs")
})
