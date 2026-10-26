#' Trait correlation at which one factor stops being an equal composite
#'
#' Two traits correlate \eqn{\varphi}, each is measured by `m` items with
#' standardized loading `lambda`, and a one-factor model is fitted to all
#' `2m` items. Above the threshold returned here the population one-factor
#' solution is an equal-weight composite of the two traits. Below it the
#' solution tilts toward one trait, and at \eqn{\varphi = 0} it represents
#' that trait alone.
#'
#' * `"ML"`: \eqn{\varphi^* = (1 + (m - 1)\lambda^2) / (m(2 - \lambda^2))}
#'   (Theorem 1 of the paper).
#' * `"ULS"`: \eqn{1/(2m)}, whatever the loading (Theorem 2).
#' * `"WLSMV"`: the threshold of diagonally weighted least squares on
#'   polychoric correlations with the thresholds free (the weighted
#'   least-squares result of Section 5.1 of the paper). It solves
#'   \eqn{\varphi = 1/2 - (m - 1)/(2m\omega(\varphi))}, where \eqn{\omega} is
#'   the ratio of the asymptotic variance of a within-block polychoric
#'   correlation to that of a cross-block one, computed from the influence
#'   function of the two-step polychoric estimator with numerically
#'   differentiated cell probabilities. The root is searched on
#'   \eqn{(0, 1/2)}; a single root was found in every design evaluated, but
#'   uniqueness is not established in general. `NA` is returned when the
#'   composite condition already holds at \eqn{\varphi = 0}.
#'
#' @param m Number of items per trait (at least 2). Vectorized.
#' @param lambda Standardized loading, common to all items. Vectorized.
#' @param estimator One of `"ML"`, `"ULS"` or `"WLSMV"`.
#' @param thresholds For `"WLSMV"` only: the item thresholds on the latent
#'   response scale, common to all items, for example `c(-0.8, 0, 0.8)` for
#'   four ordered categories.
#' @return A numeric vector of thresholds.
#' @examples
#' phi_star(m = c(2, 4, 8), lambda = 0.7)          # .493 .409 .367
#' phi_star(m = 4, lambda = 0.7, estimator = "ULS")  # .125
#' phi_star(m = 4, lambda = 0.5, estimator = "WLSMV",
#'          thresholds = c(-0.8, 0, 0.8))          # .081
#' @export
phi_star <- function(m, lambda, estimator = c("ML", "ULS", "WLSMV"),
                     thresholds = NULL) {
  estimator <- match.arg(estimator)
  check_m(m, "m")
  check_loading(lambda, "lambda")
  n <- max(length(m), length(lambda))
  m <- rep_len(m, n); lambda <- rep_len(lambda, n)
  L <- lambda^2
  switch(estimator,
    ML = (1 + (m - 1) * L) / (m * (2 - L)),
    ULS = 1 / (2 * m),
    WLSMV = {
      if (is.list(thresholds))
        stop("phi_star() assumes thresholds shared by all items; for item-specific thresholds use pop_fit(R, \"WLSMV\", thresholds = <list>).",
             call. = FALSE)
      if (is.null(thresholds) || !is.numeric(thresholds) || anyNA(thresholds) ||
          anyDuplicated(thresholds))
        stop("`thresholds` (a numeric vector of distinct values shared by all items) is required for WLSMV.",
             call. = FALSE)
      thr <- sort(thresholds)
      mapply(wlsmv_threshold, m, lambda, MoreArgs = list(thr = thr))
    })
}

kappa_wlsmv <- function(L, phi, thr) {
  avar_polychoric(L, thr, thr) / avar_polychoric(L * phi, thr, thr)
}

phi_kappa <- function(m, kappa) 0.5 - (m - 1) / (2 * m * kappa)

wlsmv_threshold <- function(m, lambda, thr) {
  L <- lambda^2
  g <- function(phi) phi - phi_kappa(m, kappa_wlsmv(L, phi, thr))
  if (g(0) >= 0) return(NA_real_)
  stats::uniroot(g, c(0, 0.5), tol = 1e-10)$root
}
