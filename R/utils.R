## Internal helpers

check_loading <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || any(abs(x) <= 0 | abs(x) >= 1))
    stop(sprintf("`%s` must contain loadings with 0 < |loading| < 1.", name), call. = FALSE)
  invisible(x)
}

check_m <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || any(x < 2) || any(x != round(x)))
    stop(sprintf("`%s` must be a whole number of at least 2.", name), call. = FALSE)
  invisible(x)
}

check_phi <- function(x, name = "phi") {
  if (!is.numeric(x) || anyNA(x) || any(abs(x) >= 1))
    stop(sprintf("`%s` must lie strictly between -1 and 1.", name), call. = FALSE)
  invisible(x)
}

## Run `expr` with a fixed seed and restore the caller's random number state.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  expr
}

## Root mean square residual over the p(p + 1)/2 nonredundant elements
srmr_of <- function(R, Sigma) {
  res <- (R - Sigma)[lower.tri(R, diag = TRUE)]
  sqrt(mean(res^2))
}

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
