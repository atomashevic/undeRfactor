test_that("the quartic agrees with numerical minimization", {
  set.seed(3)
  for (i in 1:25) {
    mA <- sample(2:9, 1); mB <- sample(2:9, 1)
    lA <- runif(1, 0.3, 0.9); lB <- runif(1, 0.3, 0.9); phi <- runif(1, 0.01, 0.95)
    b <- block_solutions(mA, mB, lA, lB, phi)
    f <- pop_fit(two_trait_cor(phi, lA, lB, mA, mB), exact = FALSE, nstart = 10)
    expect_equal(min(b$F), f$F, tolerance = 1e-8)
    expect_equal(sum(b$type == "minimum"), nrow(f$minima))
  }
})

test_that("folds match Proposition 2 and the threshold matches Theorem 1", {
  expect_equal(fold_phi(4, 4, 0.68, 0.72), 0.2416, tolerance = 1e-3)
  expect_equal(fold_phi(4, 4, 0.60, 0.80), 0.0165, tolerance = 1e-2)
  expect_equal(fold_phi(4, 4, 0.7, 0.7), phi_star(4, 0.7), tolerance = 1e-6)
  ## the condition of Proposition 2 fails: no weaker-block minimum
  expect_length(fold_phi(4, 4, 0.5, 0.8), 0)
})

test_that("the exchangeable design has mirror minima below phi* and one above", {
  below <- block_solutions(4, 4, 0.7, 0.7, 0.3)
  expect_equal(sum(below$type == "minimum"), 2)
  expect_equal(sum(below$global), 2)
  expect_equal(below$F[1], 0.940795, tolerance = 1e-6)
  above <- block_solutions(4, 4, 0.7, 0.7, 0.5)
  expect_equal(nrow(above), 1)
  expect_equal(above$share_A, 0.5, tolerance = 1e-8)
})

test_that("unequal block lengths: smooth switch and discontinuous jump (Section 5.2)", {
  ## two items at .85 and three at .75: the shorter block is favored above .758
  expect_equal(block_solutions(2, 3, 0.85, 0.75, 0.70)$follows[1], "B")
  expect_equal(block_solutions(2, 3, 0.85, 0.75, 0.80)$follows[1], "A")
  ## two items at .90 and four at .70: two minima exchange global status near .326
  lo <- block_solutions(2, 4, 0.90, 0.70, 0.32)
  hi <- block_solutions(2, 4, 0.90, 0.70, 0.33)
  expect_equal(sum(lo$type == "minimum"), 2)
  expect_lt(lo$share_A[1], 0.2)
  expect_gt(hi$share_A[1], 0.9)
})

test_that("negative phi reflects block B", {
  pos <- block_solutions(9, 7, 0.7, 0.8, 0.6)
  neg <- block_solutions(9, 7, 0.7, 0.8, -0.6)
  expect_equal(neg$F, pos$F)
  expect_equal(neg$magnitude_share_A, pos$share_A)
  expect_equal(neg$l_B, -pos$l_B)
})

test_that("block_solutions checks its input", {
  expect_error(block_solutions(4, 4, 0.7, 0.7, 0), "pop_fit")
  expect_error(block_solutions(4, 4, c(0.6, 0.7), 0.7, 0.3), "single values")
})

test_that("negative generating loadings keep the sign pattern of the population", {
  for (d in list(c(0.3, -0.7, 0.6), c(-0.3, 0.7, -0.6), c(-0.4, -0.7, -0.6))) {
    b <- block_solutions(4, 3, d[2], d[3], d[1])
    f <- pop_fit(two_trait_cor(d[1], d[2], d[3], m_A = 4, m_B = 3), exact = FALSE, nstart = 5)
    expect_equal(min(b$F), f$F, tolerance = 1e-8)
    expect_equal(sign(b$l_A[1] * b$l_B[1]), sign(prod(d)))
    expect_equal(b$magnitude_share_A[1], f$share$magnitude_share_A, tolerance = 1e-6)
  }
})
