test_that("closed forms of Theorem 1 agree with numerical minimization", {
  for (phi in c(0.05, 0.3, 0.6, 0.9)) {
    R <- two_trait_cor(phi, 0.7, m_A = 4)
    a <- pop_fit(R)
    b <- pop_fit(R, exact = FALSE)
    expect_equal(a$F, b$F, tolerance = 1e-8)
    expect_equal(a$share$share, b$share$share, tolerance = 1e-6)
    expect_equal(nrow(a$minima), nrow(b$minima))
  }
})

test_that("population CFI is 1/2 at phi = 0", {
  for (m in c(2, 4, 8)) expect_equal(pop_fit(two_trait_cor(0, 0.6, m_A = m))$CFI, 0.5)
})

test_that("below phi* the factor follows one trait, above it is an equal composite", {
  below <- pop_fit(two_trait_cor(0.3, 0.7, m_A = 4))
  above <- pop_fit(two_trait_cor(0.5, 0.7, m_A = 4))
  expect_gt(below$share$share, 0.75)
  expect_equal(nrow(below$minima), 2)
  expect_equal(above$share$share, 0.5)
  expect_equal(above$share$follows, "neither")
})

test_that("SRMR matches its closed form above phi*", {
  f <- pop_fit(two_trait_cor(0.6, 0.7, m_A = 4))
  expect_equal(f$SRMR, 0.49 * 0.4 * sqrt(12 / 63), tolerance = 1e-10)
})

test_that("ULS matches Theorem 2 on both sides of 1/p", {
  m <- 4; L <- 0.49
  f1 <- pop_fit(two_trait_cor(0.05, 0.7, m_A = m), "ULS")
  expect_equal(f1$F, m * (m - 1) * L^2 / 2 - m^2 * (m - 1) * L^2 * 0.05^2, tolerance = 1e-10)
  f2 <- pop_fit(two_trait_cor(0.3, 0.7, m_A = m), "ULS")
  expect_equal(f2$F, m^2 * (m - 1) * L^2 * 0.7^2 / (2 * m - 1), tolerance = 1e-10)
  expect_equal(f2$share$share, 0.5, tolerance = 1e-6)
})

test_that("WLSMV tilts below its threshold and not above it", {
  thr <- c(-0.8, 0, 0.8)
  expect_equal(round(pop_fit(two_trait_cor(0.06, 0.5, m_A = 4), "WLSMV", thr)$share$share, 3), 0.715)
  expect_equal(pop_fit(two_trait_cor(0.10, 0.5, m_A = 4), "WLSMV", thr)$share$share, 0.5,
               tolerance = 1e-5)
})

test_that("Proposition 1: the weaker block keeps a local minimum below the fold", {
  low <- pop_fit(two_trait_cor(0.1, lambda_A = 0.68, lambda_B = 0.72))
  expect_equal(nrow(low$minima), 2)
  expect_equal(low$share$follows, "B")
  expect_gt(low$minima$share_A[2], 0.5)
  high <- pop_fit(two_trait_cor(0.3, lambda_A = 0.68, lambda_B = 0.72))
  expect_equal(nrow(high$minima), 1)
})

test_that("narrow basins near the fold are found (.60/.80, fold at .017)", {
  near <- pop_fit(two_trait_cor(0.01, lambda_A = 0.60, lambda_B = 0.80))
  expect_equal(nrow(near$minima), 2)
  expect_equal(round(near$minima$F, 4), c(0.6064, 1.9922))
  past <- pop_fit(two_trait_cor(0.03, lambda_A = 0.60, lambda_B = 0.80))
  expect_equal(nrow(past$minima), 1)
})

test_that("PHQ-ADS shares of Section 7 are reproduced", {
  share <- function(phi, lam, est = "ML", thr = NULL)
    pop_fit(two_trait_cor(phi, lam, m_A = 9, m_B = 7), est, thresholds = thr)$share$share_A
  expect_equal(round(c(share(0.653, 0.70), share(0.653, 0.80)), 2), c(0.64, 0.70))
  expect_equal(round(c(share(0.912, 0.70), share(0.912, 0.80)), 2), c(0.58, 0.59))
  ## the depression share depends on the estimator: .665 under ML, .568 under WLSMV
  expect_equal(share(0.653, 0.75), 0.665, tolerance = 1e-3)
  expect_equal(share(0.653, 0.75, "WLSMV", c(-0.8, 0, 0.8)), 0.568, tolerance = 1e-3)
  esc <- pop_fit(two_trait_cor(0.865, 0.75, m_A = 9, m_B = 7))
  expect_equal(round(c(esc$CFI, esc$RMSEA), 3), c(0.953, 0.065))
})

test_that("bad input is rejected", {
  expect_error(pop_fit(diag(3)), "at least 4")
  expect_error(pop_fit(matrix(0.5, 4, 4)), "unit diagonal")
})

test_that("the optimizer survives starts that previously left the bounds", {
  ## the PHQ-ADS layout at loading .70 and phi = .8 aborted in 0.2.0 drafts
  R <- two_trait_cor(0.8, 0.7, 0.7, m_A = 9, m_B = 7)
  exact <- pop_fit(R)
  num <- pop_fit(R, exact = FALSE)
  expect_equal(exact$method, "exact (block quartic)")
  expect_equal(num$F, exact$F, tolerance = 1e-8)
  expect_equal(exact$share$share_A, 0.602, tolerance = 1e-3)
})

test_that("negative phi: minima are sign invariant and bipolar shares are reported", {
  bip <- pop_fit(two_trait_cor(-0.6, 0.7, m_A = 4))
  pos <- pop_fit(two_trait_cor(0.6, 0.7, m_A = 4))
  expect_equal(nrow(bip$minima), 1)
  expect_equal(bip$F, pos$F)
  expect_equal(bip$minima$magnitude_share_A, 0.5, tolerance = 1e-8)
  expect_lt(bip$share$w_A * bip$share$w_B, 0)
  expect_output(print(bip), "bipolar contrast")
  expect_equal(nrow(pop_fit(two_trait_cor(-0.6, 0.7, m_A = 4), exact = FALSE)$minima), 1)
  tilt <- pop_fit(two_trait_cor(-0.3, 0.7, m_A = 4), exact = FALSE)
  expect_equal(nrow(tilt$minima), 2)
  expect_equal(tilt$minima$magnitude_share_A[1], 0.786, tolerance = 1e-3)
  expect_output(print(tilt), "mirror images")
  het <- function(phi) pop_fit(two_trait_cor(phi, c(0.6, 0.7, 0.8, 0.75), c(0.55, 0.65, 0.7, 0.8)))
  expect_equal(het(-0.5)$F, het(0.5)$F, tolerance = 1e-8)
  expect_equal(het(-0.5)$share$magnitude_share_A, het(0.5)$share$share_A, tolerance = 1e-6)
})

test_that("score shares of Sections 5.2 and 7 are reproduced", {
  f <- pop_fit(two_trait_cor(0.3, 0.7, m_A = 4))
  expect_equal(round(f$scores$share_A, 3), c(0.750, 0.500))
  sc <- function(phi, lam) pop_fit(two_trait_cor(phi, lam, m_A = 9, m_B = 7))$scores$share_A[1]
  expect_equal(round(c(sc(0.912, 0.70), sc(0.912, 0.80)), 2), c(0.57, 0.58))
  expect_equal(round(c(sc(0.653, 0.70), sc(0.653, 0.80)), 2), c(0.63, 0.68))
})
