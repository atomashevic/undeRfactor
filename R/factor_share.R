#' Covariance-matching coefficients and shares of two traits
#'
#' Solves the block-average covariance equations
#' \deqn{l_A = \lambda_A (w_A + \varphi w_B), \qquad l_B = \lambda_B (\varphi w_A + w_B).}
#' With `mode = "signed"`, the signs supplied by the caller are retained.
#' The default `mode = "absolute"` preserves the original API: it replaces
#' all four loadings by their absolute values before solving. Those weights
#' match the absolute-value targets, not generally the signed input targets.
#' A call with negative loadings and no explicit `mode` warns about this choice.
#' Trait orientations and `phi` must be consistent; reflecting one trait
#' requires reflecting its reference loading and `phi` together.
#' Signed shares depend on the chosen trait orientations; magnitude shares
#' are invariant to reflecting either trait coordinate.
#'
#' `share_A` is the signed coefficient ratio \eqn{w_A/(w_A+w_B)}. It is
#' invariant to a global sign reversal of the fitted factor, but can be outside
#' \eqn{[0, 1]} when the coefficients have opposite signs. `magnitude_share_A` is
#' \eqn{|w_A|/(|w_A|+|w_B|)}, which lies in \eqn{[0, 1]}. Neither quantity is a
#' variance share or the composition of a regression, Bartlett or sum score.
#' The covariance-matching composite need not have unit variance. With
#' heterogeneous items these equations match block averages only.
#'
#' Nonzero reference loadings and \eqn{|\varphi|<1} are required to identify
#' both coefficients. A zero reference loading removes an equation (or is
#' incompatible with a nonzero fitted loading). At \eqn{|\varphi|=1}, the
#' trait covariance matrix is singular; the two coefficients are not unique
#' even if the target is compatible. Such inputs are rejected.
#'
#' If both fitted loadings are zero, all shares and `follows` are `NA` and
#' `status` is `"zero_composite"`. If the target ratios sum to exactly zero
#' but are nonzero, the signed shares and `follows` are `NA`, the magnitude
#' shares remain defined, and `status` is `"zero_sum"`. No tolerance truncates
#' small nonzero sums or correlations near the boundary.
#'
#' The calculation uses sum and contrast coordinates and computes `share_A`
#' directly, without summing large cancelling coefficients. This avoids
#' unnecessary roundoff but cannot remove sensitivity inherent in the inputs.
#' `condition_number` is the spectral condition number of the trait covariance
#' matrix, \eqn{(1+|\varphi|)/(1-|\varphi|)}. `sum_condition` is
#' \eqn{(|a|+|b|)/|a+b|}, where \eqn{a=l_A/\lambda_A} and
#' \eqn{b=l_B/\lambda_B} after the chosen preprocessing; it diagnoses
#' cancellation in the signed-share denominator. It is infinite for a nonzero
#' zero-sum target and `NA` for the zero composite. These are diagnostics,
#' not automatic rejection thresholds or statistical uncertainty estimates.
#'
#' @param l_A,l_B Average standardized one-factor loadings of blocks A and B.
#'   Zero is allowed. All numerical arguments must be nonempty finite real
#'   numeric vectors, either scalar or of one common length. Scalars are
#'   recycled; other length mismatches are rejected.
#' @param lambda_A,lambda_B Average standardized two-factor reference loadings,
#'   each nonzero. Inputs are used as supplied: loading magnitudes are not
#'   clipped or screened for model admissibility.
#' @param phi Trait correlation, strictly between -1 and 1.
#' @param mode `"absolute"` (the backward-compatible default) or `"signed"`.
#' @return A data frame whose first five columns retain the original API:
#'   `w_A`, `w_B`, `share_A`, `share` (the larger signed coefficient ratio,
#'   possibly greater than one), and `follows` (`"A"`, `"B"`, or `"neither"`
#'   when `abs(share_A - .5) < 1e-6`, the original descriptive tie tolerance).
#'   `follows` compares signed ratios, not coefficient signs or variances.
#'   Additional columns are `magnitude_share_A`, `magnitude_share` (the larger
#'   magnitude share), `condition_number`, `sum_condition`, `status`
#'   (`"ok"`, `"zero_sum"`, or `"zero_composite"`), and `mode`.
#' @examples
#' fit <- pop_fit(two_trait_cor(phi = 0.653, lambda_A = 0.75, m_A = 9, m_B = 7))
#' fit$share
#' factor_share(0.70, 0.62, 0.75, 0.75, 0.65)
#' factor_share(0.6, -0.2, 0.7, 0.7, 0.3, mode = "signed")
#' factor_share(0.5, -0.5, 0.7, 0.7, 0.3, mode = "signed")
#' @export
factor_share <- function(l_A, l_B, lambda_A, lambda_B, phi,
                         mode = c("absolute", "signed")) {
  explicit_mode <- !missing(mode)
  mode <- match.arg(mode)
  inputs <- list(l_A = l_A, l_B = l_B, lambda_A = lambda_A,
                 lambda_B = lambda_B, phi = phi)
  for (nm in names(inputs)) {
    x <- inputs[[nm]]
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
        !length(x) || any(!is.finite(x)))
      stop(sprintf("`%s` must be a nonempty finite real numeric vector.", nm),
           call. = FALSE)
  }
  n <- max(lengths(inputs))
  if (any(!lengths(inputs) %in% c(1L, n)))
    stop("Arguments must be scalar or have one common length.", call. = FALSE)
  check_phi(phi)
  if (any(lambda_A == 0) || any(lambda_B == 0))
    stop("Reference loadings must be nonzero: a zero reference makes the coefficient system unidentified or incompatible.",
         call. = FALSE)
  inputs <- lapply(inputs, rep_len, length.out = n)
  if (mode == "absolute") {
    if (!explicit_mode && any(unlist(inputs[1:4], use.names = FALSE) < 0))
      warning('Default mode = "absolute" discards loading signs; use mode = "signed" for signed covariance matching, or specify mode = "absolute" explicitly.',
              call. = FALSE)
    inputs[1:4] <- lapply(inputs[1:4], abs)
  }
  a <- inputs$l_A / inputs$lambda_A
  b <- inputs$l_B / inputs$lambda_B
  if (any(!is.finite(a)) || any(!is.finite(b)))
    stop("Loading ratios overflow floating-point range.", call. = FALSE)
  phi <- inputs$phi
  ## Scale before sums, differences and absolute-value normalization.
  scale <- pmax(abs(a), abs(b))
  zero <- scale == 0
  scale[zero] <- 1
  a <- a / scale; b <- b / scale
  target_sum <- a + b
  sum_part <- (target_sum / 2) / (1 + phi)
  contrast_part <- ((a - b) / 2) / (1 - phi)
  u <- sum_part + contrast_part
  v <- sum_part - contrast_part
  w_A <- scale * u; w_B <- scale * v
  if (any(!is.finite(w_A)) || any(!is.finite(w_B)))
    stop("Covariance-matching coefficients overflow floating-point range.", call. = FALSE)
  defined <- target_sum != 0
  share_A <- rep(NA_real_, n)
  share_A[defined] <- 0.5 + ((a[defined] - b[defined]) / target_sum[defined]) *
    ((1 + phi[defined]) / (2 * (1 - phi[defined])))
  if (any(!is.finite(share_A[defined])))
    stop("Signed coefficient shares overflow floating-point range.", call. = FALSE)
  follows <- ifelse(is.na(share_A), NA_character_,
                    ifelse(abs(share_A - 0.5) < 1e-6, "neither",
                           ifelse(share_A > 0.5, "A", "B")))
  magnitude_share_A <- rep(NA_real_, n)
  magnitude_share_A[!zero] <- abs(u[!zero]) / (abs(u[!zero]) + abs(v[!zero]))
  sum_condition <- (abs(a) + abs(b)) / abs(target_sum)
  sum_condition[zero] <- NA_real_
  data.frame(w_A = w_A, w_B = w_B, share_A = share_A,
             share = pmax(share_A, 1 - share_A), follows = follows,
             magnitude_share_A = magnitude_share_A,
             magnitude_share = pmax(magnitude_share_A, 1 - magnitude_share_A),
             condition_number = (1 + abs(phi)) / (1 - abs(phi)),
             sum_condition = sum_condition,
             status = ifelse(zero, "zero_composite", ifelse(defined, "ok", "zero_sum")),
             mode = mode, stringsAsFactors = FALSE)
}

## Fit callers retain usable estimates even when a share cannot be computed.
unavailable_factor_share <- function(reason) {
  data.frame(w_A = NA_real_, w_B = NA_real_, share_A = NA_real_, share = NA_real_,
             follows = NA_character_, magnitude_share_A = NA_real_,
             magnitude_share = NA_real_, condition_number = NA_real_,
             sum_condition = NA_real_, status = "unavailable", mode = "signed",
             reason = reason, stringsAsFactors = FALSE)
}

available_factor_share <- function(l_A, l_B, lambda_A, lambda_B, phi) {
  tryCatch({
    out <- factor_share(l_A, l_B, lambda_A, lambda_B, phi, mode = "signed")
    out$reason <- NA_character_
    out
  }, error = function(e) unavailable_factor_share(conditionMessage(e)))
}
