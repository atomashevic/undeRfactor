test_that("Table 1 of the paper is reproduced", {
  tab <- rbind(blind_spot(0.7, m_A = 2), blind_spot(0.7, m_A = 4), blind_spot(0.7, m_A = 8))
  expect_equal(round(tab[, "CFI"], 2), c(0.78, 0.82, 0.84))
  expect_equal(round(tab[, "TLI"], 2), c(0.87, 0.85, 0.86))
  expect_equal(round(tab[, "RMSEA"], 2), c(0.92, 0.88, 0.84))
  expect_equal(round(tab[, "SRMR"], 2), c(0.55, 0.63, 0.65))
  small <- blind_spot(0.5, m_A = 2)
  expect_equal(unname(small["SRMR"]), 0)        # passes at every phi
  expect_equal(round(unname(small["CFI"]), 2), 0.71)
})

test_that("at the returned correlation the index sits on its cutoff", {
  b <- blind_spot(0.6, m_A = 3, index = c("CFI", "RMSEA"))
  expect_equal(pop_fit(two_trait_cor(b[["CFI"]], 0.6, m_A = 3))$CFI, 0.95, tolerance = 1e-6)
  expect_equal(pop_fit(two_trait_cor(b[["RMSEA"]], 0.6, m_A = 3))$RMSEA, 0.06, tolerance = 1e-6)
})

test_that("unequal designs run through the numerical path", {
  b <- blind_spot(0.68, 0.72, m_A = 4, index = "CFI")
  expect_true(b > 0.78 && b < 0.85)
})
