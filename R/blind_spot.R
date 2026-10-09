#' Trait correlations at which fit-index cutoffs stop rejecting one factor
#'
#' For a two-trait design, finds the trait correlation above which the
#' population value of each fit index passes its cutoff when a one-factor
#' model is fitted by maximum likelihood. Above this correlation the cutoff
#' does not reject the one-factor model even with unlimited data: the
#' interval from the returned value to 1 is the cutoff's blind spot. A value
#' of 0 means that the cutoff passes the one-factor model at every
#' nonnegative trait correlation.
#'
#' The population indices are those of [pop_fit()]. The search evaluates a
#' grid of trait correlations and refines the first crossing with
#' [stats::uniroot()], so it assumes that each index improves as the traits
#' become more correlated, which holds in every design we have examined.
#'
#' @inheritParams two_trait_cor
#' @param index Fit indices to evaluate: any of `"CFI"`, `"TLI"`, `"RMSEA"`
#'   and `"SRMR"`.
#' @param cutoff Cutoffs, in the order of `index`. Defaults to the Hu and
#'   Bentler (1999) values: CFI and TLI \eqn{\ge .95}, RMSEA \eqn{\le .06},
#'   SRMR \eqn{\le .08}.
#' @param step Grid step for the trait correlation.
#' @param nstart Random starts per fit for designs without closed forms.
#' @return A named numeric vector with one trait correlation per index.
#' @references Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit
#'   indexes in covariance structure analysis: Conventional criteria versus
#'   new alternatives. *Structural Equation Modeling, 6*(1), 1-55.
#' @examples
#' ## Table 1 of the paper, loading .70 and 8 items (4 per trait)
#' blind_spot(lambda_A = 0.7, m_A = 4)
#'
#' ## a stricter CFI cutoff
#' blind_spot(lambda_A = 0.7, m_A = 4, index = "CFI", cutoff = 0.97)
#' @export
blind_spot <- function(lambda_A, lambda_B = lambda_A, m_A = 4, m_B = m_A,
                       index = c("CFI", "TLI", "RMSEA", "SRMR"), cutoff = NULL,
                       step = 0.01, nstart = 5) {
  index <- match.arg(index, c("CFI", "TLI", "RMSEA", "SRMR"), several.ok = TRUE)
  default_cut <- c(CFI = 0.95, TLI = 0.95, RMSEA = 0.06, SRMR = 0.08)
  if (is.null(cutoff)) cutoff <- default_cut[index]
  if (length(cutoff) != length(index)) stop("`cutoff` must have one value per index.", call. = FALSE)
  names(cutoff) <- index
  higher_better <- index %in% c("CFI", "TLI")
  ## force evaluation of the defaults before lambda_A and m_A are changed
  lambda_B <- lambda_B; m_B <- m_B
  values <- function(phi) {
    f <- pop_fit(two_trait_cor(phi, lambda_A, lambda_B, m_A, m_B), nstart = nstart)
    unlist(f[index])
  }
  margin <- function(v) ifelse(higher_better, v - cutoff, cutoff - v)  # >= 0 passes
  grid <- seq(0, 1 - step, by = step)
  M <- vapply(grid, function(phi) margin(values(phi)), numeric(length(index)))
  M <- matrix(M, nrow = length(index))
  out <- vapply(seq_along(index), function(k) {
    pass <- which(M[k, ] >= 0)
    if (!length(pass)) return(NA_real_)
    first <- pass[1]
    if (first == 1) return(0)
    g <- function(phi) margin(values(phi))[k]
    stats::uniroot(g, c(grid[first - 1], grid[first]), tol = 1e-8)$root
  }, 0)
  names(out) <- index
  out
}
