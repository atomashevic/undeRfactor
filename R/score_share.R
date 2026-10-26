#' Trait composition of factor scores and sum scores
#'
#' A score \eqn{S = g^\top x} computed from the items has trait part
#' \eqn{(\Lambda^\top g)^\top \theta}, where \eqn{\Lambda} holds the two-factor
#' loadings and \eqn{\theta} the two traits (Section 5.2 of the paper). The
#' share of trait A is \eqn{t_A/(t_A + t_B)} with \eqn{t = \Lambda^\top g}.
#' Regression and Bartlett scores of the one-factor model have weights
#' proportional to \eqn{\Psi^{-1} l} and therefore the same share and the
#' same correlations with the traits; the unit-weighted sum score has
#' weights equal to the item standard deviations in the standardized metric.
#'
#' Unlike the covariance-matching share of [factor_share()], the score share
#' is defined item by item, so it applies to heterogeneous loadings without
#' block averaging. The correlations of each score with the two traits are
#' computed in the two-factor population with unit item variances,
#' \eqn{\Lambda\Phi\Lambda^\top + \Psi_0}. When the traits are highly
#' correlated, both correlations can be close to one even if the share
#' departs from .5.
#'
#' Signed shares depend on the orientation of the traits. With a negative
#' trait correlation they can lie outside zero to one; `magnitude_share_A`
#' then equals the share obtained after reflecting trait B (negating the
#' second column of `Lambda` and `phi`).
#'
#' @param loadings,uniquenesses Standardized loadings and uniquenesses of the
#'   one-factor solution. Uniquenesses must be positive for the factor score.
#' @param Lambda Standardized two-factor loadings, a matrix with one row per
#'   item and one column per trait.
#' @param phi Trait correlation, strictly between -1 and 1.
#' @param item_sd Item standard deviations, to describe the sum of the raw
#'   item scores. The default, unit standard deviations, describes the sum of
#'   the standardized items.
#' @return A data frame with rows `"factor"` (regression or Bartlett score)
#'   and `"sum"` (unit-weighted sum score) and columns `score`, `share_A`
#'   (signed ratio of the trait coefficients), `magnitude_share_A`, `cor_A`
#'   and `cor_B` (correlations of the score with each trait).
#' @examples
#' ## four items per trait, loading .7, phi = .3: the factor score's share of
#' ## the followed trait is .750, against .500 for the sum score (Section 5.2)
#' f <- pop_fit(two_trait_cor(0.3, 0.7, m_A = 4))
#' f$scores
#' @export
score_share <- function(loadings, uniquenesses, Lambda, phi, item_sd = NULL) {
  Lambda <- as.matrix(Lambda)
  p <- length(loadings)
  if (!is.numeric(Lambda) || ncol(Lambda) != 2 || nrow(Lambda) != p || anyNA(Lambda))
    stop("`Lambda` must be a numeric matrix with one row per item and two columns.",
         call. = FALSE)
  if (length(uniquenesses) != p || anyNA(c(loadings, uniquenesses)))
    stop("`loadings` and `uniquenesses` must have one value per item.", call. = FALSE)
  check_phi(phi)
  if (is.null(item_sd)) item_sd <- rep(1, p)
  if (length(item_sd) != p || any(!is.finite(item_sd)) || any(item_sd <= 0))
    stop("`item_sd` must hold one positive value per item.", call. = FALSE)
  Phi <- matrix(c(1, phi, phi, 1), 2)
  Sigma <- Lambda %*% Phi %*% t(Lambda)
  diag(Sigma) <- 1
  one <- function(label, g) {
    if (is.null(g)) return(data.frame(score = label, share_A = NA_real_,
                                      magnitude_share_A = NA_real_, cor_A = NA_real_,
                                      cor_B = NA_real_, stringsAsFactors = FALSE))
    tt <- as.vector(crossprod(Lambda, g))
    if (sum(tt) < 0) { tt <- -tt; g <- -g }
    cv <- as.vector(Phi %*% tt) / sqrt(as.numeric(crossprod(g, Sigma %*% g)))
    data.frame(score = label, share_A = if (sum(tt) == 0) NA_real_ else tt[1] / sum(tt),
               magnitude_share_A = abs(tt[1]) / sum(abs(tt)), cor_A = cv[1], cor_B = cv[2],
               stringsAsFactors = FALSE)
  }
  g_factor <- if (all(uniquenesses > 0)) loadings / uniquenesses else NULL
  rbind(one("factor", g_factor), one("sum", item_sd))
}
