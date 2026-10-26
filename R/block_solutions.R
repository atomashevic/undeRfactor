#' Exact one-factor solutions of a two-block design
#'
#' For two blocks with one loading per block (`m_A` items loading `lambda_A`
#' on trait A, `m_B` items loading `lambda_B` on trait B, trait correlation
#' `phi`), every local minimum of the one-factor ML discrepancy has constant
#' fitted loadings within each block (Lemma 1 of the paper). Its positive
#' block-uniform stationary points are the positive roots of the quartic of
#' Appendix D. `block_solutions()` returns all of them, classified as minima
#' or saddles, so the global minimum and any competing local minimum are
#' exact rather than best found. `fold_phi()` returns the trait correlations
#' at which the number of minima changes.
#'
#' A negative `phi` is handled by reflecting block B (Section 5.3): the
#' solutions are those of `|phi|` with the block-B loadings negated. Signed
#' shares keep the generating trait orientation, as in [pop_fit()], and can
#' lie outside zero to one; the magnitude shares equal the shares of the
#' reflected design. With equal blocks, the change in the number of minima
#' is the symmetric transition of Theorem 1 at `phi_star()`, not a fold.
#' Within about \eqn{10^{-9}} of that threshold or of a fold, the
#' classification of nearly coincident roots is limited by floating-point
#' precision.
#'
#' The roots and their classification are exact up to polynomial root-finding
#' error; Lemma 1 makes the lowest minimum the global minimum over the
#' unrestricted one-factor model. Zero-loading and sign-changing saddles of the
#' full likelihood are not enumerated.
#'
#' @inheritParams two_trait_cor
#' @param lambda_A,lambda_B Standardized loadings, one value per block.
#' @param phi Trait correlation, nonzero and strictly between -1 and 1. At
#'   zero the solutions retain one block; use [pop_fit()].
#' @return For `block_solutions()`, a data frame with one row per stationary
#'   point, minima first and lowest discrepancy first: `type` (`"minimum"`
#'   or `"saddle"`), `global`, the ML discrepancy `F`, `CFI`, `TLI`, `RMSEA`,
#'   `SRMR`, fitted block loadings `l_A` and `l_B`, uniquenesses `psi_A` and
#'   `psi_B`, the signed coefficient share `share_A` and `magnitude_share_A`
#'   (see [factor_share()]), `follows`, and the profile coordinate `v` of
#'   Appendix D. For
#'   `fold_phi()`, the sorted trait correlations in (0, 1) at which the number
#'   of minima changes; empty if it never does.
#' @examples
#' ## Proposition 2: the weaker block keeps a local minimum below the fold
#' block_solutions(m_A = 4, m_B = 4, lambda_A = 0.68, lambda_B = 0.72, phi = 0.1)
#' fold_phi(4, 4, 0.68, 0.72)          # .242
#'
#' ## the PHQ-ADS layout: 9 and 7 items
#' block_solutions(9, 7, 0.70, 0.80, phi = 0.865)
#' @export
block_solutions <- function(m_A, m_B = m_A, lambda_A, lambda_B = lambda_A, phi) {
  check_m(m_A, "m_A"); check_m(m_B, "m_B")
  check_loading(lambda_A, "lambda_A"); check_loading(lambda_B, "lambda_B")
  check_phi(phi)
  if (length(c(m_A, m_B, lambda_A, lambda_B, phi)) != 5)
    stop("All arguments must be single values.", call. = FALSE)
  if (phi == 0)
    stop("At phi = 0 the solutions retain one block; use pop_fit().", call. = FALSE)
  ## Reflecting items changes only signs: solve the design with positive
  ## loadings and correlation, then restore sign(l_A l_B) = sign(lambda_A lambda_B phi).
  st <- quartic_stationary(m_A, m_B, abs(lambda_A), abs(lambda_B), abs(phi))
  sA <- sign(lambda_A); sB <- sign(lambda_B) * sign(phi)
  R <- unclass_cor(two_trait_cor(phi, lambda_A, lambda_B, m_A, m_B))
  F_B <- -as.numeric(determinant(R, logarithm = TRUE)$modulus)
  p <- m_A + m_B; df <- p * (p - 3) / 2; df_B <- p * (p - 1) / 2
  out <- do.call(rbind, lapply(seq_len(nrow(st)), function(i) {
    lA <- sA * st$alpha[i]; lB <- sB * st$beta[i]
    l <- c(rep(lA, m_A), rep(lB, m_B))
    Sigma <- tcrossprod(l); diag(Sigma) <- 1
    sh <- available_factor_share(lA, lB, lambda_A, lambda_B, phi)
    Fm <- max(st$F[i], 0)
    data.frame(type = st$type[i], F = st$F[i], CFI = 1 - Fm / F_B,
               TLI = 1 - (Fm / df) / (F_B / df_B), RMSEA = sqrt(Fm / df),
               SRMR = srmr_of(R, Sigma), l_A = lA, l_B = lB,
               psi_A = 1 - st$alpha[i]^2, psi_B = 1 - st$beta[i]^2,
               share_A = sh$share_A, magnitude_share_A = sh$magnitude_share_A,
               follows = sh$follows, v = st$v[i], stringsAsFactors = FALSE)
  }))
  out <- out[order(out$type != "minimum", round(out$F, 12)), ]
  out$global <- out$type == "minimum" & abs(out$F - min(out$F[out$type == "minimum"])) < 1e-10
  rownames(out) <- NULL
  out[, c("type", "global", setdiff(names(out), c("type", "global")))]
}

#' @rdname block_solutions
#' @param n Number of grid points used to locate changes before refining
#'   them by bisection.
#' @export
fold_phi <- function(m_A, m_B = m_A, lambda_A, lambda_B = lambda_A, n = 400) {
  check_m(m_A, "m_A"); check_m(m_B, "m_B")
  check_loading(lambda_A, "lambda_A"); check_loading(lambda_B, "lambda_B")
  nmin <- function(phi)
    sum(quartic_stationary(m_A, m_B, abs(lambda_A), abs(lambda_B), phi)$type == "minimum")
  ## a geometric grid near zero catches folds that occur at very small phi
  grid <- sort(unique(c(10^seq(-6, -2, length.out = 25), seq(0.01, 0.9999, length.out = n))))
  k <- vapply(grid, nmin, 0)
  ch <- which(diff(k) != 0)
  vapply(ch, function(i) {
    a <- grid[i]; b <- grid[i + 1]; ka <- k[i]
    for (it in 1:60) { mid <- (a + b) / 2; if (nmin(mid) == ka) a <- mid else b <- mid }
    (a + b) / 2
  }, 0)
}

## Positive roots of the Appendix D quartic and their profile values; phi > 0.
quartic_stationary <- function(mA, mB, lamA, lamB, phi) {
  LA <- lamA^2; LB <- lamB^2; aA <- 1 - LA; aB <- 1 - LB
  bA <- 1 + (mA - 1) * LA; bB <- 1 + (mB - 1) * LB; nu <- mA + mB - 1
  c <- lamA * lamB * phi; kA <- (mA - 1) * aA; kB <- (mB - 1) * aB; cA <- mA * c; cB <- mB * c
  e0 <- mA * bB + mB * bA - bA * bB - cA * cB
  J <- nu * bA - e0; K <- nu * bB - e0; f <- nu - kB; g <- nu - kA
  P <- c(bA * cA * kA, cA^2 * g - bA * J + cA * cB * kA,
         bB * cA * g + cA * K - bA * cB * f - cB * J,
         bB * K - cA * cB * kB - cB^2 * f, -bB * cB * kB)
  ## Near multiple roots polyroot() returns real roots with small imaginary
  ## parts (about 1e-6 next to the triple root at phi*). Accept a candidate
  ## when its real part leaves a negligible residual.
  r <- polyroot(P)
  pw <- function(x) x^(seq_along(P) - 1)
  v <- Re(r[Re(r) > 0 & abs(Im(r)) <= 1e-4 * Mod(r)])
  v <- v[vapply(v, function(x) abs(sum(P * pw(x))) <= 1e-8 * sum(abs(P) * pw(x)), TRUE)]
  v <- sort(v)
  v <- v[c(TRUE, diff(v) > 1e-6 * v[-1])]
  G <- function(v) {
    U <- bA + cB * v; V <- bB + cA / v; H <- e0 + cB * kB * v + cA * kA / v
    nu * log(H) - (mB - 1) * log(U) - (mA - 1) * log(V) - nu * log(nu) -
      log(bA * bB - mA * mB * c^2) - (mA - 1) * log(aA) - (mB - 1) * log(aB)
  }
  ## G'(v) has the sign of -P(v), so a root is a minimum when P decreases
  ## through it. A double root is not a minimum; a triple root (the
  ## exchangeable threshold) is, with zero curvature (Appendix D).
  deriv <- function(P, x, k) {
    j <- seq_along(P) - 1
    co <- P * vapply(j, function(e) if (e < k) 0 else prod(e - seq_len(k) + 1), 0)
    sum(co * x^pmax(j - k, 0))
  }
  out <- do.call(rbind, lapply(v, function(x) {
    U <- bA + cB * x; V <- bB + cA / x; H <- e0 + cB * kB * x + cA * kA / x
    scale <- sum(abs(P) * pw(x))
    d1 <- deriv(P, x, 1) * x / scale; d2 <- deriv(P, x, 2) * x^2 / scale
    type <- if (abs(d1) > 1e-8) { if (d1 < 0) "minimum" else "saddle" } else
      if (abs(d2) > 1e-6) "saddle" else if (deriv(P, x, 3) < 0) "minimum" else "saddle"
    data.frame(v = x, F = G(x), alpha = sqrt(1 - H / (nu * V)), beta = sqrt(1 - H / (nu * U)),
               type = type, stringsAsFactors = FALSE)
  }))
  ## Within about 1e-9 of a degenerate point the classification is limited by
  ## floating-point precision; the minimum discrepancy always exists.
  if (!any(out$type == "minimum")) out$type[which.min(out$F)] <- "minimum"
  out
}
