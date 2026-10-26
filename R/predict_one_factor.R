#' Predict the one-factor solution from a two-factor fit
#'
#' From standardized two-factor estimates (loadings and trait correlation),
#' computes the population one-factor solution at four levels of detail.
#' Comparing them with the one-factor model fitted to the same data shows
#' how much of the solution a few summary numbers determine.
#'
#' * `T0`: within-block total correlations only. At \eqn{\varphi = 0} the
#'   factor retains the block with the larger total correlation and
#'   CFI = max(TC_A, TC_B)/(TC_A + TC_B) (Section 5.2; for loadings that vary
#'   within blocks this rule is a numerical finding, not a theorem).
#' * `T1`: the closed forms of Theorem 1 with one pooled loading. Defined only
#'   for equal block sizes.
#' * `T2`: one summary loading per block, solved exactly by the quartic of
#'   Appendix D ([block_solutions()]); `fold` gives the trait correlations at
#'   which the number of minima changes.
#' * `T3`: a numerical fit to the full two-factor-implied correlation matrix.
#'   It uses no theorem and serves as the reference for loading
#'   heterogeneity: if the two-factor model fits, T3 and the observed
#'   one-factor fit agree almost by construction.
#'
#' `bound_F` is the upper bound of Appendix E on the ML discrepancy, valid for
#' heterogeneous loadings.
#'
#' Shares follow the conventions of [pop_fit()]: with a negative trait
#' correlation the signed share can lie outside zero to one, and
#' `magnitude_share_A` equals the share after reflecting block B.
#'
#' @param x A fitted lavaan two-factor model with simple structure, or a
#'   numeric matrix of standardized loadings with one row per item and two
#'   columns (zero for items outside a block).
#' @param phi Trait correlation; taken from `x` when it is a lavaan fit.
#' @param estimator `"ML"`, `"ULS"` or `"WLSMV"` (see [pop_fit()]).
#' @param thresholds For `"WLSMV"`: item thresholds, as in [pop_fit()].
#' @param lambda_summary How a block's loadings are summarized for T1 and T2:
#'   `"mean"`, `"snr"` (the loading whose signal-to-noise ratio
#'   \eqn{\lambda^2/(1-\lambda^2)} is the block average), or `"tc"` (the
#'   common loading with the same total correlation as the block).
#' @param nstart,seed Passed to [pop_fit()] for T3.
#' @return An object of class `"one_factor_prediction"`: a list with `tiers`
#'   (one row per tier with the discrepancy, fit indices, block-mean fitted
#'   loadings, coefficient share, factor-score share and number of minima),
#'   `blocks` (items, summary loading and total correlation of each block),
#'   `phi`, `T0`, `fold` and `bound_F`.
#' @examples
#' Lambda <- cbind(c(0.80, 0.70, 0.60, 0.65, 0, 0, 0, 0),
#'                 c(0, 0, 0, 0, 0.75, 0.70, 0.65, 0.60))
#' predict_one_factor(Lambda, phi = 0.3)
#' @export
predict_one_factor <- function(x, phi = NULL, estimator = c("ML", "ULS", "WLSMV"),
                               thresholds = NULL, lambda_summary = c("mean", "snr", "tc"),
                               nstart = 20, seed = 1) {
  estimator <- match.arg(estimator)
  lambda_summary <- match.arg(lambda_summary)
  S_within <- NULL
  if (inherits(x, "lavaan")) {
    st <- lavaan::lavInspect(x, "std")
    if (ncol(st$lambda) != 2) stop("The lavaan model must have two factors.", call. = FALSE)
    Lambda <- st$lambda
    phi <- st$psi[1, 2]
    ss <- lavaan::lavInspect(x, "sampstat")$cov
    S_within <- stats::cov2cor(ss[rownames(Lambda), rownames(Lambda)])
  } else {
    Lambda <- as.matrix(x)
    if (is.null(phi)) stop("`phi` is required when `x` is a loading matrix.", call. = FALSE)
  }
  check_phi(phi)
  if (!is.numeric(Lambda) || ncol(Lambda) != 2 || anyNA(Lambda))
    stop("`x` must give a numeric loading matrix with two columns.", call. = FALSE)
  inA <- Lambda[, 1] != 0; inB <- Lambda[, 2] != 0
  if (any(inA & inB) || any(!inA & !inB))
    stop("Each item must load on exactly one factor (simple structure).", call. = FALSE)
  lamA <- Lambda[inA, 1]; lamB <- Lambda[inB, 2]
  if (any(c(lamA, lamB) < 0))
    stop("Negative loadings: reverse-score those items before predicting.", call. = FALSE)
  if (any(c(lamA, lamB) >= 1))
    stop("A standardized loading is at least 1 (Heywood case); predictions need |loading| < 1.",
         call. = FALSE)
  mA <- length(lamA); mB <- length(lamB)
  check_m(mA, "items in block A"); check_m(mB, "items in block B")

  tc_of <- function(R) -0.5 * as.numeric(determinant(R, logarithm = TRUE)$modulus)
  block_R <- function(l) { R <- tcrossprod(l); diag(R) <- 1; R }
  TC <- if (is.null(S_within)) c(tc_of(block_R(lamA)), tc_of(block_R(lamB))) else
    c(tc_of(S_within[inA, inA]), tc_of(S_within[inB, inB]))
  summarize <- function(l) switch(lambda_summary,
    mean = mean(l),
    snr = { s <- mean(l^2 / (1 - l^2)); sqrt(s / (1 + s)) },
    tc = {
      m <- length(l); target <- tc_of(block_R(l))
      stats::uniroot(function(z) -0.5 * ((m - 1) * log(1 - z^2) + log(1 + (m - 1) * z^2)) -
                       target, c(1e-8, 1 - 1e-8), tol = 1e-12)$root
    })
  lbar <- c(summarize(lamA), summarize(lamB))

  thr_for <- function(idx) if (is.list(thresholds)) thresholds[idx] else thresholds
  row_of <- function(tier, f) {
    sc <- if (is.null(f$scores)) NA_real_ else f$scores$share_A[1]
    data.frame(tier = tier, method = f$method, F = f$F, CFI = f$CFI, TLI = f$TLI,
               RMSEA = f$RMSEA, SRMR = f$SRMR, l_A = f$minima$l_A[1], l_B = f$minima$l_B[1],
               share_A = f$share$share_A, magnitude_share_A = f$share$magnitude_share_A,
               follows = f$share$follows, score_share_A = sc,
               n_minima = nrow(f$minima), stringsAsFactors = FALSE)
  }
  fit_tier <- function(R, thr) pop_fit(R, estimator, thresholds = thr, nstart = nstart, seed = seed)
  ordA <- c(which(inA), which(inB))
  tiers <- list()
  if (mA == mB) {
    pooled <- summarize(c(lamA, lamB))
    tiers$T1 <- row_of("T1", fit_tier(two_trait_cor(phi, pooled, m_A = mA), thr_for(ordA)))
  }
  tiers$T2 <- row_of("T2", fit_tier(two_trait_cor(phi, lbar[1], lbar[2], mA, mB), thr_for(ordA)))
  tiers$T3 <- row_of("T3", fit_tier(two_trait_cor(phi, lamA, lamB), thr_for(ordA)))
  tiers <- do.call(rbind, tiers)
  rownames(tiers) <- NULL

  ## Appendix E: bound on F from the two-factor parameters
  Phi <- matrix(c(1, phi, phi, 1), 2)
  Lm <- matrix(0, mA + mB, 2); Lm[seq_len(mA), 1] <- lamA; Lm[mA + seq_len(mB), 2] <- lamB
  e <- eigen(Phi, symmetric = TRUE)
  Ph <- e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% t(e$vectors)
  kap <- sort(eigen(Ph %*% crossprod(Lm, Lm / (1 - rowSums(Lm^2))) %*% Ph, symmetric = TRUE,
                    only.values = TRUE)$values, decreasing = TRUE)
  p <- mA + mB
  bound_F <- (p - 1) * log(1 + kap[2] / (p - 1)) - log(1 + kap[2])

  blocks <- data.frame(block = c("A", "B"), items = c(mA, mB), lambda_bar = lbar, TC = TC,
                       stringsAsFactors = FALSE)
  structure(list(tiers = tiers, blocks = blocks, phi = phi,
                 T0 = list(follows_at_zero = c("A", "B")[which.max(TC)],
                           CFI_at_zero = max(TC) / sum(TC)),
                 fold = if (estimator == "ML") fold_phi(mA, mB, lbar[1], lbar[2]) else NULL,
                 bound_F = if (estimator == "ML") bound_F else NA_real_,
                 estimator = estimator, lambda_summary = lambda_summary,
                 items = list(A = rownames(Lambda)[inA], B = rownames(Lambda)[inB])),
            class = "one_factor_prediction")
}

#' @export
print.one_factor_prediction <- function(x, ...) {
  cat(sprintf("Predicted one-factor solution (%s), %d + %d items, phi = %s\n", x$estimator,
              x$blocks$items[1], x$blocks$items[2], fmt(x$phi)))
  cat(sprintf("  Block summaries (%s): %s (A), %s (B); total correlations %s and %s.\n",
              x$lambda_summary, fmt(x$blocks$lambda_bar[1]), fmt(x$blocks$lambda_bar[2]),
              fmt(x$blocks$TC[1]), fmt(x$blocks$TC[2])))
  show <- x$tiers[, c("tier", "F", "CFI", "RMSEA", "SRMR", "l_A", "l_B", "share_A",
                      "magnitude_share_A", "score_share_A", "n_minima")]
  num <- vapply(show, is.numeric, TRUE)
  show[num] <- lapply(show[num], function(v) round(v, 4))
  print(show, row.names = FALSE)
  if (length(x$fold))
    cat(sprintf("  T2: the number of minima changes at phi = %s.\n",
                paste(fmt(x$fold), collapse = ", ")))
  if (is.finite(x$bound_F))
    cat(sprintf("  Upper bound on F for these loadings (Appendix E): %s.\n", fmt(x$bound_F, 4)))
  invisible(x)
}
