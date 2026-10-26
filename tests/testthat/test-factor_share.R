test_that("factor_share inverts the composite loadings", {
  lam_A <- 0.7; lam_B <- 0.6; phi <- 0.4; w <- c(0.8, 0.3)
  l_A <- lam_A * (w[1] + phi * w[2]); l_B <- lam_B * (phi * w[1] + w[2])
  s <- factor_share(l_A, l_B, lam_A, lam_B, phi)
  expect_equal(c(s$w_A, s$w_B), w)
  expect_equal(s$share_A, 0.8 / 1.1)
  expect_equal(s$follows, "A")
})

test_that("an equal composite has share .5 and signs do not matter", {
  s <- factor_share(-0.5, 0.5, 0.7, 0.7, 0.3, mode = "absolute")
  expect_equal(s$share, 0.5)
  expect_equal(s$follows, "neither")
})

test_that("signed coefficients reconstruct signed block covariance targets", {
  cases <- expand.grid(l_A = c(-.6, 0, .6), l_B = c(-.2, 0, .2),
                       lambda_A = c(-.7, .7), lambda_B = c(-.6, .6),
                       phi = c(-.8, 0, .8))
  s <- do.call(factor_share, c(cases, list(mode = "signed")))
  expect_equal(cases$lambda_A * (s$w_A + cases$phi * s$w_B), cases$l_A)
  expect_equal(cases$lambda_B * (cases$phi * s$w_A + s$w_B), cases$l_B)
  reflected <- cases
  reflected$l_A <- -cases$l_A; reflected$l_B <- -cases$l_B
  r <- do.call(factor_share, c(reflected, list(mode = "signed")))
  expect_equal(r$w_A, -s$w_A); expect_equal(r$w_B, -s$w_B)
  expect_equal(r$share_A, s$share_A)
  expect_equal(r$magnitude_share_A, s$magnitude_share_A)
  ## Consistent change of trait A coordinates (not a fitted-factor reflection).
  reflected <- cases
  reflected$lambda_A <- -cases$lambda_A; reflected$phi <- -cases$phi
  r <- do.call(factor_share, c(reflected, list(mode = "signed")))
  expect_equal(r$w_A, -s$w_A); expect_equal(r$w_B, s$w_B)
  expect_equal(r$magnitude_share_A, s$magnitude_share_A)
})

test_that("legacy preprocessing is explicit and does not reconstruct opposite signs", {
  expect_warning(s <- factor_share(.6, -.2, .7, .7, .3), "discards loading signs")
  expect_no_warning(a <- factor_share(.6, -.2, .7, .7, .3, mode = "absolute"))
  expect_equal(s, a)
  z <- factor_share(.6, -.2, .7, .7, .3, mode = "signed")
  expect_equal(.7 * (.3 * a$w_A + a$w_B), .2)
  expect_equal(.7 * (.3 * z$w_A + z$w_B), -.2)
  expect_gt(z$share_A, 1)
  expect_lt(z$magnitude_share_A, 1)
  ## Positive targets can also imply opposite-sign coefficients.
  z <- factor_share(.7, .07, .7, .7, .5)
  expect_lt(z$w_B, 0)
  expect_gt(z$share_A, 1)
  expect_equal(z$magnitude_share_A, abs(z$w_A) / (abs(z$w_A) + abs(z$w_B)))
})

test_that("zero coefficient sums do not masquerade as equal composites", {
  s <- factor_share(c(0, .5, 0), c(0, -.5, .5), .7, .7, .3, mode = "signed")
  expect_equal(s$status, c("zero_composite", "zero_sum", "ok"))
  expect_true(all(is.na(s$share_A[1:2])))
  expect_true(all(is.na(s$share[1:2])))
  expect_true(all(is.na(s$follows[1:2])))
  expect_true(is.na(s$magnitude_share_A[1]))
  expect_equal(s$magnitude_share_A[2], .5)
  expect_equal(s$w_A[1], 0); expect_equal(s$w_B[1], 0)
  expect_true(is.na(s$sum_condition[1]))
  expect_equal(s$sum_condition[2], Inf)
  expect_error(factor_share(.5, .5, 0, .7, .3), "unidentified or incompatible")
  expect_error(factor_share(0, .5, 0, .7, .3), "unidentified or incompatible")
})

test_that("near-boundary results retain intrinsic contrast amplification", {
  phi <- c(.9, .99, .9999, 1 - 1e-8, 1 - .Machine$double.eps / 2,
           -.9, -1 + .Machine$double.eps / 2)
  equal <- factor_share(.6, .6, 1, 1, phi, mode = "signed")
  expect_equal(equal$share_A, rep(.5, length(phi)))
  expect_equal(equal$w_A, .6 / (1 + phi))
  expect_equal(equal$condition_number, (1 + abs(phi)) / (1 - abs(phi)))
  a <- .6 + 1e-6; b <- .6
  z <- factor_share(a, b, 1, 1, phi, mode = "signed")
  ## Direct identity is independent of the eigen-coordinate implementation.
  expected <- (a - phi * b) / ((1 - phi) * (a + b))
  expect_equal(z$share_A, expected, tolerance = 1e-9)
  expect_gt(z$share_A[5], 1e8)
  ## At this boundary, summing huge opposite weights can lose the denominator.
  expect_true(is.finite(z$share_A[5]))
  near_cancel <- factor_share(.5, -.5 + 2^-50, 1, 1, .3, mode = "signed")
  expect_equal(near_cancel$status, "ok")
  expect_gt(near_cancel$sum_condition, 1e14)
  expect_gt(near_cancel$share_A, 1e14)
  for (p in c(-1, 1, -1.01, 1.01))
    expect_error(factor_share(.5, .5, .7, .7, p), "strictly between")
})

test_that("validation rejects malformed inputs without partial recycling", {
  for (nm in c("l_A", "l_B", "lambda_A", "lambda_B", "phi")) {
    for (bad in list(numeric(), NA_real_, NaN, Inf, -Inf, "0.5", TRUE,
                     1 + 1i, matrix(.5))) {
      args <- list(l_A = .5, l_B = .4, lambda_A = .7, lambda_B = .7, phi = .3)
      args[[nm]] <- bad
      expect_error(do.call(factor_share, args), "finite real numeric vector")
    }
  }
  expect_error(factor_share(c(.5, .6), c(.1, .2, .3), .7, .7, .3), "common length")
  expect_error(factor_share(.5, .4, .7, .7, .3, mode = "wrong"), "arg")
  expect_equal(nrow(factor_share(c(.5, .6), .4, .7, .7, .3)), 2)
  expect_equal(factor_share(2, 2, 1, 1, 0)$w_A, 2)
  expect_error(factor_share(1, 1, 1e-320, 1, .3), "ratios overflow")
})

test_that("well-conditioned legacy results and scale invariance are retained", {
  cases <- expand.grid(l_A = c(0, .3, .6), l_B = c(0, .2, .5),
                       lambda_A = c(.6, .8), lambda_B = c(.6, .8), phi = c(-.3, 0, .7))
  cases <- cases[cases$l_A + cases$l_B > 0, ]
  s <- do.call(factor_share, cases)
  a <- abs(cases$l_A / cases$lambda_A); b <- abs(cases$l_B / cases$lambda_B)
  w_A <- (a - cases$phi * b) / (1 - cases$phi^2)
  w_B <- (b - cases$phi * a) / (1 - cases$phi^2)
  expect_equal(s$w_A, w_A); expect_equal(s$w_B, w_B)
  expect_equal(s$share_A, w_A / (w_A + w_B))
  expect_identical(names(s)[1:5], c("w_A", "w_B", "share_A", "share", "follows"))
  base <- factor_share(.6, -.2, .7, .7, .3, mode = "signed")
  for (scale in c(1e-280, 1e280)) {
    z <- factor_share(scale * .6, scale * -.2, .7, .7, .3, mode = "signed")
    expect_equal(z$share_A, base$share_A)
    expect_equal(z$magnitude_share_A, base$magnitude_share_A)
    expect_equal(z$w_A / scale, base$w_A)
    expect_equal(z$w_B / scale, base$w_B)
  }
})
