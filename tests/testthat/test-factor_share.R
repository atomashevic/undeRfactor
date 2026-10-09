test_that("factor_share inverts the composite loadings", {
  lam_A <- 0.7; lam_B <- 0.6; phi <- 0.4; w <- c(0.8, 0.3)
  l_A <- lam_A * (w[1] + phi * w[2]); l_B <- lam_B * (phi * w[1] + w[2])
  s <- factor_share(l_A, l_B, lam_A, lam_B, phi)
  expect_equal(c(s$w_A, s$w_B), w)
  expect_equal(s$share_A, 0.8 / 1.1)
  expect_equal(s$follows, "A")
})

test_that("an equal composite has share .5 and signs do not matter", {
  s <- factor_share(-0.5, 0.5, 0.7, 0.7, 0.3)
  expect_equal(s$share, 0.5)
  expect_equal(s$follows, "neither")
})
