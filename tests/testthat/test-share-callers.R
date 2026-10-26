test_that("block summaries preserve both fitted and reference signs", {
  d <- attr(two_trait_cor(-.3, -.7, .6, m_A = 2, m_B = 2), "two_trait")
  l <- c(.6, .6, -.2, -.2)
  s <- share_of(l, d)
  expect_equal(c(s$l_A, s$l_B), c(.6, -.2))
  expect_equal(s$share, available_factor_share(.6, -.2, -.7, .6, -.3))
  r <- share_of(-l, d)$share
  expect_equal(r$w_A, -s$share$w_A); expect_equal(r$share_A, s$share_A)
  zero <- share_of(rep(0, 4), d)$share
  expect_equal(zero$status, "zero_composite")
  one <- share_of(c(0, 0, .2, .2), d)$share
  expect_equal(one$status, "ok")
  expect_true(is.finite(one$w_A))
})

test_that("unavailable reference shares retain an explanatory result", {
  for (ref in list(c(0, .7, .3), c(.7, 0, .3), c(.7, .7, 1),
                   c(.7, .7, -1), c(NA, .7, .3), c(.7, .7, Inf))) {
    s <- available_factor_share(.5, .4, ref[1], ref[2], ref[3])
    expect_equal(s$status, "unavailable")
    expect_true(is.na(s$share_A))
    expect_true(nzchar(s$reason))
  }
  d <- attr(two_trait_cor(.3, c(-.7, .7), c(-.6, .6)), "two_trait")
  s <- share_of(c(.5, -.5, .4, -.4), d)
  expect_equal(s$share$status, "unavailable")
})

test_that("population fits pass their signed block targets through", {
  for (R in list(two_trait_cor(-.3, .7, .6, m_A = 3),
                 two_trait_cor(.3, -.7, .6, m_A = 3))) {
    fit <- pop_fit(R, nstart = 2)
    d <- attr(R, "two_trait")
    target <- c(mean(fit$loadings[d$blocks$A]), mean(fit$loadings[d$blocks$B]))
    expect_lt(prod(target), 0)
    expect_equal(c(fit$minima$l_A[1], fit$minima$l_B[1]), target)
    expect_equal(c(mean(d$lambda_A) * (fit$share$w_A + d$phi * fit$share$w_B),
                   mean(d$lambda_B) * (d$phi * fit$share$w_A + fit$share$w_B)), target)
    expect_equal(fit$share$mode, "signed")
  }
  fit <- pop_fit(two_trait_cor(.3, c(-.7, .7), c(-.6, .6)), nstart = 2)
  expect_true(is.finite(fit$F))
  expect_equal(fit$share$status, "unavailable")
  expect_true(all(fit$minima$share_status == "unavailable"))
  expect_output(print(fit), "Trait share unavailable")
})

test_that("order checks preserve signed solutions and survive bad references", {
  skip_on_cran()
  set.seed(20261011)
  R <- two_trait_cor(-.4, .7, .6, m_A = 3)
  X <- as.data.frame(matrix(rnorm(900 * 6), 900) %*% chol(R))
  names(X) <- rownames(R)
  blocks <- list(A = paste0("A", 1:3), B = paste0("B", 1:3))
  fit <- order_check(X, blocks, nstart = 0)
  for (i in seq_len(nrow(fit$solutions))) {
    l <- lavaan::lavInspect(fit$fits[[i]], "std")$lambda[, 1][names(X)]
    if (sum(l) < 0) l <- -l
    target <- c(mean(l[1:3]), mean(l[4:6]))
    expect_lt(prod(target), 0)
    expect_equal(unname(unlist(fit$solutions[i, c("l_A", "l_B")])), target)
    ref <- fit$two_factor
    direct <- factor_share(target[1], target[2], ref$lambda_A, ref$lambda_B,
                           ref$phi, mode = "signed")
    expect_equal(fit$solutions$share_A[i], direct$share_A)
  }
  ## Make only the two-factor reference invalid; one-factor fits still succeed.
  original <- lavaan::lavInspect
  bad_kind <- "boundary"
  testthat::local_mocked_bindings(lavInspect = function(object, what, ...) {
    st <- original(object, "std")
    two <- ncol(st$lambda) == 2
    if (two && bad_kind == "unconverged" && what == "converged") return(FALSE)
    if (two && bad_kind == "extraction" && what == "std") stop("Reference extraction failed.")
    ans <- original(object, what, ...)
    if (two && identical(what, "std")) {
      if (bad_kind == "boundary") ans$psi[1, 2] <- 1
      if (bad_kind == "zero") ans$lambda[, 1] <- 0
      if (bad_kind == "nonfinite") ans$lambda[, 1] <- NA_real_
    }
    ans
  }, .package = "lavaan")
  for (kind in c("boundary", "zero", "nonfinite", "unconverged", "extraction")) {
    bad_kind <- kind
    bad <- order_check(X, blocks, nstart = 0)
    expect_gt(bad$n_converged, 0)
    expect_true(all(is.na(bad$solutions$share_A)))
    expect_true(all(bad$solutions$share_status == "unavailable"))
    expect_true(all(nzchar(bad$solutions$share_reason)))
  }
})
