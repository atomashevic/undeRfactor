test_that("ML threshold matches Theorem 1", {
  expect_equal(phi_star(c(2, 4, 8), 0.7), c(1.49 / 3.02, 0.4089404, 0.3667219), tolerance = 1e-6)
  expect_equal(phi_star(c(4, 8), 0.5), c(0.25, 2.75 / 14))
})

test_that("ULS threshold is 1/p whatever the loading", {
  expect_equal(phi_star(c(2, 4, 8), c(0.3, 0.7, 0.9), "ULS"), c(0.25, 0.125, 0.0625))
})

test_that("WLSMV thresholds match Section 3 of the paper", {
  sym <- c(-0.8, 0, 0.8); skew <- c(0.3, 1.0, 1.7)
  expect_equal(phi_star(c(4, 8), 0.5, "WLSMV", sym), c(0.081, 0.010), tolerance = 0.001 / 0.081)
  expect_equal(round(phi_star(c(4, 8), 0.5, "WLSMV", sym), 3), c(0.081, 0.010))
  expect_equal(round(phi_star(c(4, 8), 0.5, "WLSMV", skew), 3), c(0.076, 0.004))
  expect_equal(round(phi_star(2, 0.7, "WLSMV", sym), 3), 0.098)
  expect_true(all(is.na(phi_star(c(4, 8), 0.7, "WLSMV", sym))))
  expect_error(phi_star(4, 0.5, "WLSMV"), "thresholds")
})

test_that("inputs are checked", {
  expect_error(phi_star(1, 0.5), "at least 2")
  expect_error(phi_star(4, 1.2), "loadings")
})

test_that("item-specific WLSMV thresholds are refused with a pointer to pop_fit", {
  expect_error(phi_star(4, 0.7, "WLSMV", thresholds = list(c(-1, 0, 1), c(0, 1))), "pop_fit")
})
