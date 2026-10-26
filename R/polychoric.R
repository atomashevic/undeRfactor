## Asymptotic variance of the two-step polychoric correlation estimator
## (thresholds from the margins, then rho by bivariate ML given the
## thresholds), from its influence function with numerically differentiated
## cell probabilities. Used for the WLSMV weights.

bvn_cdf <- function(x, y, rho) {
  out <- numeric(length(x))
  fin <- is.finite(x) & is.finite(y)
  out[x == -Inf | y == -Inf] <- 0
  both_inf <- x == Inf & y == Inf
  out[both_inf] <- 1
  xi <- x == Inf & is.finite(y); out[xi] <- stats::pnorm(y[xi])
  yi <- y == Inf & is.finite(x); out[yi] <- stats::pnorm(x[yi])
  if (any(fin)) {
    out[fin] <- if (abs(rho) < 1e-14) stats::pnorm(x[fin]) * stats::pnorm(y[fin]) else
      pbivnorm::pbivnorm(x[fin], y[fin], rho = rho)
  }
  out
}

## Cell probabilities of two ordinal items with thresholds t1, t2
cell_probs <- function(rho, t1, t2) {
  a <- c(-Inf, t1, Inf); b <- c(-Inf, t2, Inf)
  g <- expand.grid(i = seq_along(a), j = seq_along(b))
  Fm <- matrix(bvn_cdf(a[g$i], b[g$j], rho), length(a), length(b))
  Fm[-1, -1] - Fm[-length(a), -1] - Fm[-1, -length(b)] + Fm[-length(a), -length(b)]
}

## N times the asymptotic variance of the polychoric estimate of rho
avar_polychoric <- function(rho, t1, t2, h = 1e-5) {
  P <- cell_probs(rho, t1, t2)
  d_rho <- (cell_probs(rho + h, t1, t2) - cell_probs(rho - h, t1, t2)) / (2 * h)
  s_rho <- d_rho / P
  I_rr <- sum(P * s_rho^2)
  K1 <- length(t1); K2 <- length(t2)
  corr <- 0 * P
  for (a in seq_len(K1)) {
    e <- replace(0 * t1, a, h)
    s_t <- ((cell_probs(rho, t1 + e, t2) - cell_probs(rho, t1 - e, t2)) / (2 * h)) / P
    IF_t <- outer((seq_len(K1 + 1) <= a) - stats::pnorm(t1[a]), rep(1, K2 + 1)) /
      stats::dnorm(t1[a])
    corr <- corr + sum(P * s_rho * s_t) * IF_t
  }
  for (b in seq_len(K2)) {
    e <- replace(0 * t2, b, h)
    s_t <- ((cell_probs(rho, t1, t2 + e) - cell_probs(rho, t1, t2 - e)) / (2 * h)) / P
    IF_t <- outer(rep(1, K1 + 1), (seq_len(K2 + 1) <= b) - stats::pnorm(t2[b])) /
      stats::dnorm(t2[b])
    corr <- corr + sum(P * s_rho * s_t) * IF_t
  }
  IF <- (s_rho - corr) / I_rr
  sum(P * IF^2)
}

## Weight matrix (inverse asymptotic variances) for the off-diagonal elements
## of a polychoric correlation matrix R. `thresholds` is one vector shared by
## all items or a list with one vector per item.
dwls_weights <- function(R, thresholds) {
  p <- nrow(R)
  thr <- if (is.list(thresholds)) thresholds else rep(list(thresholds), p)
  if (length(thr) != p) stop("`thresholds` must be one vector or a list with one vector per item.",
                             call. = FALSE)
  key_thr <- vapply(thr, function(t) paste(signif(t, 12), collapse = ","), "")
  cache <- new.env(hash = TRUE)
  W <- matrix(0, p, p)
  for (i in 2:p) for (j in seq_len(i - 1)) {
    key <- paste(signif(R[i, j], 12), key_thr[i], key_thr[j], sep = "|")
    key2 <- paste(signif(R[i, j], 12), key_thr[j], key_thr[i], sep = "|")
    w <- cache[[key]]
    if (is.null(w)) w <- cache[[key2]]
    if (is.null(w)) {
      w <- 1 / avar_polychoric(R[i, j], thr[[i]], thr[[j]])
      cache[[key]] <- w
    }
    W[i, j] <- W[j, i] <- w
  }
  W
}
