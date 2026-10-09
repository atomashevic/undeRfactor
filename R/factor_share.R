#' Share of each trait in a fitted single factor
#'
#' Writes the fitted factor as a composite \eqn{w_A \theta_A + w_B \theta_B}
#' of the two traits whose covariances with the items reproduce the
#' block-average fitted loadings,
#' \deqn{l_A = \lambda_A (w_A + \varphi w_B), \qquad l_B = \lambda_B (\varphi w_A + w_B),}
#' and reports the share of trait A, \eqn{w_A/(w_A + w_B)}. The share is .5
#' for an equal composite and 1 (or 0) when the factor represents trait A
#' (or B) alone.
#'
#' The fitted loadings `l_A`, `l_B` come from the one-factor model; the true
#' or estimated loadings `lambda_A`, `lambda_B` and trait correlation `phi`
#' come from the two-factor model. All loadings are standardized, and block
#' averages are taken in absolute value, which removes sign reflections.
#'
#' @param l_A,l_B Average standardized one-factor loadings of the items of
#'   block A and block B.
#' @param lambda_A,lambda_B Average standardized loadings of the two blocks
#'   in the two-factor model.
#' @param phi Trait correlation in the two-factor model.
#' @return A data frame with the weights `w_A`, `w_B`, the share of trait A
#'   `share_A`, the share of the trait the factor follows `share`
#'   (\eqn{\ge .5}), and `follows` (`"A"`, `"B"`, or `"neither"` for an equal
#'   composite).
#' @examples
#' ## population one-factor loadings for 9 depression and 7 anxiety items
#' fit <- pop_fit(two_trait_cor(phi = 0.653, lambda_A = 0.75, m_A = 9, m_B = 7))
#' fit$share
#' factor_share(l_A = 0.70, l_B = 0.62, lambda_A = 0.75, lambda_B = 0.75, phi = 0.65)
#' @export
factor_share <- function(l_A, l_B, lambda_A, lambda_B, phi) {
  check_phi(phi)
  a <- abs(l_A) / abs(lambda_A)
  b <- abs(l_B) / abs(lambda_B)
  w_A <- (a - phi * b) / (1 - phi^2)
  w_B <- (b - phi * a) / (1 - phi^2)
  share_A <- w_A / (w_A + w_B)
  follows <- ifelse(abs(share_A - 0.5) < 1e-6, "neither",
                    ifelse(share_A > 0.5, "A", "B"))
  data.frame(w_A = w_A, w_B = w_B, share_A = share_A,
             share = pmax(share_A, 1 - share_A), follows = follows,
             stringsAsFactors = FALSE)
}
