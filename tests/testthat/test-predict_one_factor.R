Lambda <- cbind(c(0.80, 0.70, 0.60, 0.65, 0, 0, 0, 0),
                c(0, 0, 0, 0, 0.75, 0.70, 0.65, 0.60))

test_that("the tiers reproduce pop_fit and block_solutions", {
  pr <- predict_one_factor(Lambda, phi = 0.3)
  expect_equal(pr$tiers$tier, c("T1", "T2", "T3"))
  t3 <- pop_fit(two_trait_cor(0.3, Lambda[1:4, 1], Lambda[5:8, 2]))
  expect_equal(pr$tiers$F[3], t3$F)
  b <- block_solutions(4, 4, mean(Lambda[1:4, 1]), mean(Lambda[5:8, 2]), 0.3)
  expect_equal(pr$tiers$F[2], min(b$F))
  expect_equal(pr$tiers$n_minima[2], sum(b$type == "minimum"))
  expect_equal(pr$fold, fold_phi(4, 4, mean(Lambda[1:4, 1]), mean(Lambda[5:8, 2])))
  ## Appendix E: the bound holds for heterogeneous loadings
  expect_lte(pr$tiers$F[3], pr$bound_F)
})

test_that("unequal block sizes have no T1 tier", {
  L9 <- cbind(c(rep(0.7, 9), rep(0, 7)), c(rep(0, 9), rep(0.8, 7)))
  pr <- predict_one_factor(L9, phi = 0.8)
  expect_equal(pr$tiers$tier, c("T2", "T3"))
  expect_equal(pr$tiers$F[1], pr$tiers$F[2], tolerance = 1e-8)
})

test_that("a negative trait correlation gives the reflected solutions", {
  pos <- predict_one_factor(Lambda, phi = 0.3)
  neg <- predict_one_factor(Lambda, phi = -0.3)
  expect_equal(neg$tiers$F, pos$tiers$F, tolerance = 1e-8)
  expect_equal(neg$tiers$magnitude_share_A, pos$tiers$share_A, tolerance = 1e-6)
  expect_true(all(neg$tiers$l_B < 0))
})

test_that("a lavaan two-factor fit is accepted", {
  skip_on_cran()
  set.seed(2)
  R <- two_trait_cor(0.25, c(0.75, 0.68, 0.62, 0.55), c(0.73, 0.69, 0.63, 0.59))
  X <- as.data.frame(matrix(rnorm(3000 * 8), 3000) %*% chol(unclass(R)))
  names(X) <- rownames(R)
  f2 <- lavaan::cfa("fA =~ A1 + A2 + A3 + A4\nfB =~ B1 + B2 + B3 + B4", data = X, std.lv = TRUE)
  pr <- predict_one_factor(f2)
  expect_equal(pr$items$A, paste0("A", 1:4))
  expect_equal(pr$phi, lavaan::lavInspect(f2, "std")$psi[1, 2])
  expect_output(print(pr), "Predicted one-factor solution")
})

test_that("predict_one_factor checks its input", {
  expect_error(predict_one_factor(Lambda), "phi")
  cross <- Lambda; cross[1, 2] <- 0.2
  expect_error(predict_one_factor(cross, phi = 0.3), "simple structure")
  neg <- Lambda; neg[1, 1] <- -0.5
  expect_error(predict_one_factor(neg, phi = 0.3), "reverse-score")
  hey <- Lambda; hey[1, 1] <- 1.01
  expect_error(predict_one_factor(hey, phi = 0.3), "Heywood")
})
