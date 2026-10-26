#' Population one-factor fit to a correlation matrix
#'
#' Fits a one-factor model to a population correlation matrix, finds its
#' closed-form solution where available or searches numerically, and returns the population discrepancy, the fit indices,
#' the fitted loadings and, for matrices made by [two_trait_cor()], the share
#' of each trait in the fitted factor. With two traits the one-factor
#' discrepancy can have two minima (one per trait); all distinct minima found
#' are listed.
#'
#' * `"ML"`: the discrepancy is
#'   \eqn{F = \mathrm{tr}(\Sigma^{-1}R) - \log\det(\Sigma^{-1}R) - p}, which
#'   is twice lavaan's `fmin`. CFI = 1 - F/F_B, TLI = 1 - (F/df)/(F_B/df_B)
#'   and RMSEA = sqrt(F/df), with F_B = -log det R, df = p(p - 3)/2 and
#'   df_B = p(p - 1)/2. For the exchangeable design (equal loadings and block
#'   sizes) the closed forms of Theorem 1 are used. For other designs with one
#'   loading per block and \eqn{\varphi \ne 0}, all minima come from the exact
#'   quartic of Appendix D (see [block_solutions()]). Otherwise the
#'   uniquenesses are optimized numerically on the log scale from several
#'   deterministic starts and `nstart` random ones; these minima are the best
#'   found, not certified.
#' * `"ULS"`: the least-squares discrepancy over the off-diagonal elements,
#'   \eqn{\sum_{i<j} (R_{ij} - l_i l_j)^2}, with free uniquenesses.
#' * `"WLSMV"`: the same sum weighted by the inverse asymptotic variances of
#'   the polychoric correlations, which is the population WLSMV (DWLS)
#'   solution when `R` holds the polychoric correlations. Requires
#'   `thresholds`.
#'
#' For `"ULS"` and `"WLSMV"` only SRMR is returned among the fit indices:
#' their CFI, TLI and RMSEA depend on scaling choices (unscaled, scaled,
#' robust) that a population value cannot settle without them.
#'
#' Trait shares use signed block means and the generating trait orientation.
#' They match block-average covariances, not necessarily individual items or
#' factor-score composition. An invalid reference or unrepresentable result
#' leaves the fit intact and returns an unavailable share with a `reason`.
#' For \eqn{\varphi < 0} the signed coefficients of a solution can have
#' opposite signs, so `share_A` can lie far outside zero to one, or be huge
#' for a bipolar composite whose coefficients nearly cancel. The magnitude
#' share then describes the solution; it equals the share obtained after
#' reflecting block B (Section 5.3 of the paper), and the printout reports it.
#' Minima that differ only in the overall sign of the loadings are the same
#' solution; the reported sign makes the loading sum nonnegative.
#'
#' `scores` gives the trait composition of the regression (or Bartlett)
#' factor score and of the unit-weighted sum of the standardized items; see
#' [score_share()].
#'
#' @param R A population correlation matrix, for example from
#'   [two_trait_cor()].
#' @param estimator One of `"ML"`, `"ULS"` or `"WLSMV"`.
#' @param thresholds For `"WLSMV"`: item thresholds on the latent response
#'   scale, one vector shared by all items or a list with one vector per item.
#' @param nstart Number of random starts added to the deterministic ones.
#' @param exact Use the closed forms where they exist (ML, exchangeable
#'   design). Set to `FALSE` to force numerical minimization.
#' @param seed Seed for the random starts; the caller's random number state is
#'   restored afterwards. `NULL` uses the current state.
#' @return An object of class `"pop_fit"`: a list with `estimator`, `F` (the
#'   lowest discrepancy found), `F_B`, `df`, `CFI`, `TLI`, `RMSEA`, `SRMR`,
#'   `loadings`, `uniquenesses`, `share` (a row of [factor_share()], or
#'   `NULL`), `scores` (the result of [score_share()], or `NULL`), `minima`
#'   (one row per distinct minimum, lowest discrepancy first, with signed and
#'   magnitude shares) and `method`.
#' @examples
#' R <- two_trait_cor(phi = 0.3, lambda_A = 0.7, m_A = 4)
#' pop_fit(R)                    # below phi* = .409: follows one trait
#' pop_fit(two_trait_cor(0.6, 0.7, m_A = 4))   # above phi*: equal composite
#'
#' ## unequal loadings: the weaker block keeps a spurious local minimum
#' pop_fit(two_trait_cor(0.1, lambda_A = 0.68, lambda_B = 0.72))$minima
#' @export
pop_fit <- function(R, estimator = c("ML", "ULS", "WLSMV"), thresholds = NULL,
                    nstart = 20, exact = TRUE, seed = 1) {
  estimator <- match.arg(estimator)
  R <- check_cor(R)
  p <- nrow(R)
  design <- attr(R, "two_trait")
  W <- NULL
  if (estimator == "WLSMV") {
    if (is.null(thresholds))
      stop("`thresholds` is required for WLSMV.", call. = FALSE)
    W <- dwls_weights(unclass_cor(R), thresholds)
  }
  if (estimator == "ULS") W <- matrix(1, p, p) - diag(p)

  if (estimator == "ML" && exact && is_exchangeable(design)) {
    sols <- ml_exchangeable(design)
    method <- "closed form"
  } else if (estimator == "ML" && exact && is_block_uniform(design)) {
    sols <- ml_blocks(design)
    method <- "exact (block quartic)"
  } else {
    sols <- with_seed(seed, {
      if (estimator == "ML") ml_numeric(unclass_cor(R), design, nstart) else
        ls_numeric(unclass_cor(R), W, design, nstart)
    })
    method <- "numerical"
  }
  sols <- lapply(distinct_minima(sols), function(s) { s$l <- orient(s$l); s })
  best <- sols[[1]]
  l <- best$l; psi <- best$psi
  names(l) <- names(psi) <- rownames(R)
  Sigma <- tcrossprod(l) + diag(psi, p)

  F_B <- -as.numeric(determinant(unclass_cor(R), logarithm = TRUE)$modulus)
  df <- p * (p - 3) / 2; df_B <- p * (p - 1) / 2
  idx <- c(CFI = NA_real_, TLI = NA_real_, RMSEA = NA_real_)
  if (estimator == "ML") {
    Fm <- max(best$F, 0)
    idx <- c(CFI = 1 - Fm / F_B, TLI = 1 - (Fm / df) / (F_B / df_B), RMSEA = sqrt(Fm / df))
  }
  minima <- do.call(rbind, lapply(sols, function(s) {
    sh <- share_of(s$l, design)
    data.frame(F = s$F, l_A = sh$l_A, l_B = sh$l_B, share_A = sh$share_A,
               share_status = sh$share$status, share_reason = sh$share$reason,
               min_uniqueness = min(s$psi), magnitude_share_A = sh$magnitude_share_A)
  }))
  rownames(minima) <- NULL
  if (is.null(design)) minima <- minima[, c("F", "min_uniqueness")]
  share <- if (is.null(design)) NULL else share_of(best$l, design)$share
  scores <- if (is.null(design)) NULL else
    score_share(l, psi, design_lambda(design), design$phi)
  structure(list(estimator = estimator, F = best$F, F_B = F_B, df = df,
                 CFI = unname(idx["CFI"]), TLI = unname(idx["TLI"]),
                 RMSEA = unname(idx["RMSEA"]), SRMR = srmr_of(unclass_cor(R), Sigma),
                 loadings = l, uniquenesses = psi, share = share, scores = scores,
                 minima = minima, method = method),
            class = "pop_fit")
}

#' @export
print.pop_fit <- function(x, ...) {
  cat(sprintf("Population one-factor fit (%s, %s), %d items\n", x$estimator, x$method,
              length(x$loadings)))
  if (x$estimator == "ML") {
    cat(sprintf("  F = %s   CFI = %s   TLI = %s   RMSEA = %s   SRMR = %s\n",
                fmt(x$F, 4), fmt(x$CFI), fmt(x$TLI), fmt(x$RMSEA), fmt(x$SRMR)))
  } else {
    cat(sprintf("  discrepancy = %s   SRMR = %s\n", fmt(x$F, 5), fmt(x$SRMR)))
  }
  mirror <- nrow(x$minima) == 2 && abs(diff(x$minima$F)) < 1e-10
  if (!is.null(x$share)) {
    s <- x$share
    if (is.na(s$follows) && !bipolar(s)) {
      cat(sprintf("  Trait share unavailable (%s).\n", s$status))
      if (!is.na(s$reason)) cat(sprintf("  %s\n", s$reason))
    } else if (bipolar(s)) {
      ms <- s$magnitude_share_A
      if (abs(ms - 0.5) < 1e-6) {
        cat("  The factor is a bipolar contrast of the two traits with equal weights",
            "(magnitude share .500).\n")
      } else {
        cat(sprintf("  The factor is bipolar and weights trait %s more: magnitude share %s.\n",
                    if (ms > 0.5) "A" else "B", fmt(max(ms, 1 - ms))))
      }
      if (mirror) cat("  The two minima are mirror images, one per trait, with the same fit.\n")
    } else if (s$follows == "neither") {
      cat("  The factor is an equal composite of the two traits (share .500).\n")
    } else if (mirror) {
      cat(sprintf("  The factor follows one trait: share %s. The two minima are mirror images,\n", fmt(s$share)))
      cat("  one per trait, with the same fit.\n")
    } else {
      cat(sprintf("  The factor follows trait %s: share %s.\n", s$follows, fmt(s$share)))
    }
  }
  if (nrow(x$minima) > 1 && !mirror)
    cat(sprintf("  %d distinct minima; the lowest-discrepancy one found is shown. See $minima.\n",
                nrow(x$minima)))
  invisible(x)
}

## ------------------------------------------------------------- internals ----

check_cor <- function(R) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) stop("`R` must be a square matrix.", call. = FALSE)
  if (nrow(R) < 4) stop("`R` must have at least 4 items.", call. = FALSE)
  if (max(abs(R - t(R))) > 1e-10) stop("`R` must be symmetric.", call. = FALSE)
  if (max(abs(diag(R) - 1)) > 1e-10) stop("`R` must be a correlation matrix (unit diagonal).",
                                          call. = FALSE)
  ev <- eigen(unclass_cor(R), symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) stop("`R` must be positive definite.", call. = FALSE)
  if (is.null(rownames(R))) {
    nm <- paste0("x", seq_len(nrow(R)))
    dimnames(R) <- list(nm, nm)
  }
  R
}

unclass_cor <- function(R) { attr(R, "two_trait") <- NULL; R }

is_exchangeable <- function(d) {
  !is.null(d) && length(d$lambda_A) == length(d$lambda_B) &&
    length(unique(c(d$lambda_A, d$lambda_B))) == 1
}

is_block_uniform <- function(d) {
  !is.null(d) && d$phi != 0 && length(unique(d$lambda_A)) == 1 &&
    length(unique(d$lambda_B)) == 1
}

## Common sign: nonnegative loading sum. A bipolar solution has a sum of
## (nearly) zero; its first nonzero loading is then made positive.
orient <- function(l) {
  s <- sum(l)
  if (abs(s) <= 1e-10 * sum(abs(l))) s <- l[which.max(abs(l) > 1e-8)]
  if (s < 0) -l else l
}

share_of <- function(l, design) {
  if (is.null(design)) return(list(l_A = NA_real_, l_B = NA_real_, share_A = NA_real_,
                                 magnitude_share_A = NA_real_,
                                 share = list(status = NA_character_, reason = NA_character_)))
  lA <- mean(l[design$blocks$A]); lB <- mean(l[design$blocks$B])
  s <- available_factor_share(lA, lB, mean(design$lambda_A), mean(design$lambda_B),
                              design$phi)
  list(l_A = lA, l_B = lB, share_A = s$share_A, magnitude_share_A = s$magnitude_share_A,
       share = s)
}

## Signed coefficients of opposite sign: a bipolar factor, described by magnitudes
bipolar <- function(s) isTRUE(s$status %in% c("ok", "zero_sum") && s$w_A * s$w_B < 0)

## Two-factor loading matrix of a two_trait_cor() design
design_lambda <- function(design) {
  p <- length(design$lambda_A) + length(design$lambda_B)
  Lambda <- matrix(0, p, 2)
  Lambda[design$blocks$A, 1] <- design$lambda_A
  Lambda[design$blocks$B, 2] <- design$lambda_B
  Lambda
}

## Theorem 1: closed-form ML solutions of the exchangeable design; for
## phi < 0 those of |phi| with block B reflected back (Section 5.3)
ml_exchangeable <- function(d) {
  m <- length(d$lambda_A); lambda <- d$lambda_A[1]; phi <- abs(d$phi)
  sB <- c(rep(1, m), rep(if (d$phi < 0) -1 else 1, m))
  L <- lambda^2; a <- 1 - L; n <- m - 1; N <- 2 * m - 1; b <- 1 + n * L
  q <- m * L * phi / b; qs <- L / (2 - L)
  F0 <- -n * log(a) - log(b)
  one <- function(l) list(F = NA_real_, l = sB * l, psi = 1 - l^2, converged = TRUE)
  if (q < qs) {
    s2 <- L + (2 - L) * q^2
    sp <- sqrt(s2 + 2 * q); sm <- sqrt(max(s2 - 2 * q, 0))
    alpha <- (sp + sm) / 2; beta <- (sp - sm) / 2
    Fv <- F0 + n * log(1 - q^2)
    s1 <- one(c(rep(alpha, m), rep(beta, m))); s1$F <- Fv
    s2 <- one(c(rep(beta, m), rep(alpha, m))); s2$F <- Fv
    list(s1, s2)
  } else {
    km <- m * L * (1 - phi) / a
    s <- one(rep(sqrt(L * (n + m * phi) / N), 2 * m))
    s$F <- N * log(1 + km / N) - log(1 + km)
    list(s)
  }
}

## Appendix D: all minima of a design with one loading per block
ml_blocks <- function(d) {
  mA <- length(d$lambda_A); mB <- length(d$lambda_B)
  b <- block_solutions(mA, mB, d$lambda_A[1], d$lambda_B[1], d$phi)
  b <- b[b$type == "minimum", ]
  lapply(seq_len(nrow(b)), function(i) {
    l <- c(rep(b$l_A[i], mA), rep(b$l_B[i], mB))
    list(F = b$F[i], l = l, psi = 1 - l^2, converged = TRUE)
  })
}

## ML: concentrated discrepancy over the uniquenesses (as in stats::factanal)
ml_obj <- function(psi, R) {
  sc <- 1 / sqrt(psi)
  e <- eigen(R * outer(sc, sc), symmetric = TRUE, only.values = TRUE)$values[-1]
  sum(e - log(e) - 1)
}

ml_loadings <- function(psi, R) {
  sc <- 1 / sqrt(psi)
  E <- eigen(R * outer(sc, sc), symmetric = TRUE)
  l <- sqrt(psi) * E$vectors[, 1] * sqrt(max(E$values[1] - 1, 0))
  if (sum(l) < 0) l <- -l
  l
}

ml_grad <- function(psi, R) {
  l <- ml_loadings(psi, R)
  g <- tcrossprod(l) + diag(psi) - R
  diag(g) / psi^2
}

ml_numeric <- function(R, design, nstart) {
  p <- nrow(R)
  starts <- list(pmin(pmax((1 - 0.5 / p) / diag(solve(R)), 0.05), 1))
  if (!is.null(design)) {
    A <- seq_len(p) %in% design$blocks$A
    ## one block fitted, the other left to its uniquenesses: the phi = 0
    ## solution of each trait, which lies in the basin of that trait's minimum
    ## even when the basin is narrow (near the fold of Proposition 1)
    one_block <- function(B) {
      r <- R[B, B][lower.tri(R[B, B])]
      ifelse(B, pmin(pmax(1 - mean(r), 0.05), 0.95), 0.99)
    }
    starts <- c(starts, list(one_block(A), one_block(!A),
                             ifelse(A, 0.1, 0.9), ifelse(A, 0.9, 0.1)))
  }
  starts <- c(starts, replicate(nstart, stats::runif(p, 0.2, 1), simplify = FALSE))
  ## Optimize over log uniquenesses. On the raw scale, L-BFGS-B can step
  ## outside its bounds from starts with small uniquenesses, where the
  ## gradient is of order 1/psi^2.
  obj <- function(t, R) ml_obj(exp(t), R)
  grad <- function(t, R) { psi <- exp(t); ml_grad(psi, R) * psi }
  run <- function(s) {
    o <- stats::optim(log(s), obj, grad, R = R, method = "L-BFGS-B", lower = log(1e-4),
                      upper = 0, control = list(factr = 10, pgtol = 0, maxit = 10000))
    o$par <- exp(o$par)
    o
  }
  lapply(starts, function(s) tryCatch({
    o <- run(s)
    ## leave a saddle point if the optimizer stopped on one
    o2 <- run(pmin(pmax(o$par * (1 + stats::runif(p, -1e-3, 1e-3)), 1e-4), 1))
    if (o2$value < o$value - 1e-12) o <- o2
    list(F = o$value, l = ml_loadings(o$par, R), psi = o$par,
         converged = o$convergence == 0 && is.finite(o$value))
  }, error = function(e) NULL))
}

## Weighted least squares over the off-diagonal elements (ULS: unit weights)
ls_obj <- function(l, R, W) {
  E <- R - tcrossprod(l); diag(E) <- 0
  0.5 * sum(W * E^2)
}

ls_grad <- function(l, R, W) {
  E <- R - tcrossprod(l); diag(E) <- 0
  -2 * as.vector((W * E) %*% l)
}

ls_numeric <- function(R, W, design, nstart) {
  p <- nrow(R)
  e1 <- eigen(R, symmetric = TRUE)
  starts <- list(pmin(abs(e1$vectors[, 1]) * sqrt(max(e1$values[1] - 1, 0.1)), 0.95))
  if (!is.null(design)) {
    A <- seq_len(p) %in% design$blocks$A
    one_block <- function(B) {
      r <- R[B, B][lower.tri(R[B, B])]
      ifelse(B, sqrt(min(max(mean(r), 0.01), 0.95)), 0.01)
    }
    starts <- c(starts, list(one_block(A), one_block(!A), ifelse(A, 0.9, 0.1), ifelse(A, 0.1, 0.9)))
  }
  starts <- c(starts, replicate(nstart, stats::runif(p, 0, 1), simplify = FALSE))
  run <- function(s) stats::optim(s, ls_obj, ls_grad, R = R, W = W, method = "L-BFGS-B",
                                  lower = -1, upper = 1,
                                  control = list(factr = 10, pgtol = 0, maxit = 10000))
  lapply(starts, function(s) tryCatch({
    o <- run(s)
    o2 <- run(pmin(pmax(o$par + stats::runif(p, -1e-3, 1e-3), -1), 1))
    if (o2$value < o$value - 1e-14) o <- o2
    l <- o$par
    if (sum(l) < 0) l <- -l
    list(F = o$value, l = l, psi = 1 - l^2,
         converged = o$convergence == 0 && is.finite(o$value))
  }, error = function(e) NULL))
}

## Cluster solutions into distinct minima, global first. Loadings that differ
## only in their overall sign are the same solution.
distinct_minima <- function(sols, tol = 1e-3) {
  sols <- Filter(function(s) !is.null(s) && isTRUE(s$converged), sols)
  if (!length(sols)) stop("The optimizer did not converge from any start.", call. = FALSE)
  ## ties between mirror-image minima: the one with the larger first loading first
  sols <- sols[order(round(vapply(sols, `[[`, 0, "F"), 10),
                     -vapply(sols, function(s) abs(s$l[1]), 0))]
  same <- function(a, b) min(max(abs(a - b)), max(abs(a + b))) < tol
  out <- list()
  for (s in sols) {
    if (!any(vapply(out, function(o) same(o$l, s$l), TRUE))) out[[length(out) + 1]] <- s
  }
  out
}
