#' Does item order or the starting point change a one-factor solution?
#'
#' When the items come from two traits, the one-factor discrepancy can have
#' two minima, one following each trait, and lavaan's default starting values
#' depend on which item is listed first. `order_check()` fits the one-factor
#' model with each block listed first and from random starts, clusters the
#' converged solutions into distinct minima, and reports which trait each one
#' follows. It also fits the correlated two-factor model, whose loadings and
#' trait correlation define the share of each trait (see [factor_share()]).
#'
#' Read the result in two steps. If the two default fits reach different
#' solutions, item order decides which trait the factor follows. If a default
#' fit is not the solution with the lowest discrepancy, report the lowest one
#' and treat the other as a local minimum.
#'
#' @param data A data frame with the items.
#' @param blocks A list of two character vectors with the item names of each
#'   trait, for example `list(dep = paste0("phq", 1:9), anx = paste0("gad", 1:7))`.
#'   The names of the list label the traits (default `"A"` and `"B"`).
#' @param estimator Passed to [lavaan::cfa()], for example `"ML"` or `"WLSMV"`.
#' @param ordered Passed to [lavaan::cfa()]: `TRUE` treats all items as
#'   ordinal.
#' @param nstart Number of random starts.
#' @param tol Largest difference in standardized loadings for two fits to
#'   count as the same solution.
#' @param seed Seed for the random starts; the caller's random number state
#'   is restored afterwards. `NULL` uses the current state.
#' @param ... Further arguments passed to [lavaan::cfa()], for example
#'   `missing = "fiml"`.
#' @return An object of class `"order_check"`: a list with `solutions` (one
#'   row per distinct one-factor solution, lowest discrepancy first, with its
#'   fit, the average loadings of each block, the share of each trait and the
#'   starts that reached it), `default_first_A` and `default_first_B` (the
#'   solution reached from lavaan's default start with each block listed
#'   first), `order_dependent`, `two_factor` (block-average loadings and trait
#'   correlation of the two-factor model) and `fits` (the lavaan fit of each
#'   distinct solution).
#' @examples
#' \donttest{
#' set.seed(1)
#' R <- two_trait_cor(phi = 0.1, lambda_A = 0.68, lambda_B = 0.72)
#' X <- as.data.frame(matrix(rnorm(2000 * 8), 2000) %*% chol(unclass(R)))
#' names(X) <- rownames(R)
#' order_check(X, blocks = list(A = paste0("A", 1:4), B = paste0("B", 1:4)), nstart = 10)
#' }
#' @export
order_check <- function(data, blocks, estimator = "ML", ordered = FALSE, nstart = 20,
                        tol = 0.01, seed = 1, ...) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.list(blocks) || length(blocks) != 2)
    stop("`blocks` must be a list of two character vectors of item names.", call. = FALSE)
  labs <- names(blocks)
  if (is.null(labs) || any(!nzchar(labs))) labs <- c("A", "B")
  items <- c(blocks[[1]], blocks[[2]])
  miss <- setdiff(items, names(data))
  if (length(miss)) stop("Items not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyDuplicated(items)) stop("An item appears in both blocks.", call. = FALSE)
  X <- data[, items, drop = FALSE]
  iA <- seq_along(blocks[[1]]); iB <- length(iA) + seq_along(blocks[[2]])
  is_ord <- isTRUE(ordered) || is.character(ordered)
  ord_arg <- if (isTRUE(ordered)) items else if (is.character(ordered)) ordered else NULL

  fit1 <- function(order, l0 = NULL, psi0 = NULL) {
    nm <- items[order]
    rhs <- if (is.null(l0)) nm else sprintf("start(%.6f)*%s", l0[order], nm)
    model <- paste("f =~", paste(rhs, collapse = " + "))
    if (!is.null(psi0))
      model <- paste(c(model, sprintf("%s ~~ start(%.6f)*%s", nm, psi0[order], nm)), collapse = "\n")
    tryCatch(quiet_lavaan(lavaan::cfa(model, data = X[, nm, drop = FALSE], std.lv = TRUE,
                                      estimator = estimator, ordered = ord_arg, ...)),
             error = function(e) NULL)
  }
  sds <- if (is_ord) rep(1, length(items)) else vapply(X, stats::sd, 0, na.rm = TRUE)
  ord_AB <- c(iA, iB); ord_BA <- c(iB, iA)
  fits <- list(fit1(ord_AB), fit1(ord_BA))
  start_lab <- c("default, A first", "default, B first")
  rnd <- with_seed(seed, lapply(seq_len(nstart), function(s) {
    l0 <- stats::runif(length(items)) * sds
    psi0 <- if (is_ord) NULL else stats::runif(length(items), 0.2, 1) * sds^2
    fit1(ord_AB, l0, psi0)
  }))
  fits <- c(fits, rnd)
  start_lab <- c(start_lab, rep("random", nstart))

  ## two-factor reference
  m2 <- paste0("fA =~ ", paste(items[iA], collapse = " + "), "\n",
               "fB =~ ", paste(items[iB], collapse = " + "))
  f2 <- tryCatch(quiet_lavaan(lavaan::cfa(m2, data = X, std.lv = TRUE, estimator = estimator,
                                          ordered = ord_arg, ...)),
                 error = function(e) NULL)
  ref <- NULL
  if (!is.null(f2)) {
    st <- lavaan::lavInspect(f2, "std")
    sa <- sign(sum(st$lambda[items[iA], 1])); sb <- sign(sum(st$lambda[items[iB], 2]))
    ref <- list(lambda_A = abs(mean(st$lambda[items[iA], 1])),
                lambda_B = abs(mean(st$lambda[items[iB], 2])),
                phi = sa * sb * st$psi[1, 2])
  }

  sols <- lapply(seq_along(fits), function(k) {
    f <- fits[[k]]
    if (is.null(f) || !isTRUE(lavaan::lavInspect(f, "converged"))) return(NULL)
    l <- lavaan::lavInspect(f, "std")$lambda[, 1][items]
    if (sum(l) < 0) l <- -l
    list(k = k, fmin = unname(lavaan::fitMeasures(f, "fmin")), l = l)
  })
  ok <- which(!vapply(sols, is.null, TRUE))
  if (!length(ok)) stop("No one-factor fit converged.", call. = FALSE)
  sols <- sols[ok]
  sols <- sols[order(vapply(sols, `[[`, 0, "fmin"))]
  cl <- integer(0); reps <- list()
  sol_of <- integer(length(fits))
  for (s in sols) {
    hit <- which(vapply(reps, function(r) max(abs(r$l - s$l)) < tol, TRUE))
    if (length(hit)) sol_of[s$k] <- hit[1] else {
      reps[[length(reps) + 1]] <- s
      sol_of[s$k] <- length(reps)
    }
  }

  measures <- c("chisq", "df", "cfi", "rmsea", "srmr")
  rob <- c("cfi.robust", "rmsea.robust")
  tab <- do.call(rbind, lapply(seq_along(reps), function(j) {
    r <- reps[[j]]
    fm <- quiet_lavaan(lavaan::fitMeasures(fits[[r$k]]))
    get <- function(n) if (n %in% names(fm)) unname(fm[n]) else NA_real_
    lA <- abs(mean(r$l[iA])); lB <- abs(mean(r$l[iB]))
    sh <- if (is.null(ref)) list(share_A = NA_real_) else
      factor_share(lA, lB, ref$lambda_A, ref$lambda_B, ref$phi)
    found <- c(if (sol_of[1] == j) sprintf("default (%s first)", labs[1]),
               if (sol_of[2] == j) sprintf("default (%s first)", labs[2]),
               if (any(sol_of[-(1:2)] == j)) sprintf("%d random", sum(sol_of[-(1:2)] == j)))
    data.frame(solution = j, fmin = r$fmin, chisq = get("chisq"), df = get("df"),
               cfi = get("cfi"), rmsea = get("rmsea"), srmr = get("srmr"),
               cfi_robust = get("cfi.robust"), rmsea_robust = get("rmsea.robust"),
               l_A = lA, l_B = lB, share_A = sh$share_A,
               follows = if (is.na(sh$share_A)) NA_character_ else
                 if (abs(sh$share_A - 0.5) < 0.01) "neither" else
                   labs[ifelse(sh$share_A > 0.5, 1, 2)],
               found_by = paste(found, collapse = ", "), stringsAsFactors = FALSE)
  }))
  if (all(is.na(tab$cfi_robust))) tab$cfi_robust <- tab$rmsea_robust <- NULL
  dA <- if (sol_of[1] > 0) sol_of[1] else NA_integer_
  dB <- if (sol_of[2] > 0) sol_of[2] else NA_integer_
  structure(list(solutions = tab, default_first_A = dA, default_first_B = dB,
                 order_dependent = !is.na(dA) && !is.na(dB) && dA != dB,
                 two_factor = ref, labels = labs,
                 n_converged = length(ok), n_fits = length(fits),
                 fits = lapply(reps, function(r) fits[[r$k]])),
            class = "order_check")
}

#' @export
print.order_check <- function(x, ...) {
  s <- x$solutions
  cat(sprintf("One-factor fits: %d of %d converged, %d distinct solution%s.\n",
              x$n_converged, x$n_fits, nrow(s), if (nrow(s) == 1) "" else "s"))
  if (!is.null(x$two_factor))
    cat(sprintf("Two-factor model: mean loadings %s (%s) and %s (%s), trait correlation %s.\n",
                fmt(x$two_factor$lambda_A, 2), x$labels[1], fmt(x$two_factor$lambda_B, 2),
                x$labels[2], fmt(x$two_factor$phi, 2)))
  show <- s[, intersect(c("solution", "fmin", "chisq", "cfi", "rmsea", "share_A", "follows",
                          "found_by"), names(s))]
  names(show)[names(show) == "share_A"] <- paste0("share_", x$labels[1])
  num <- vapply(show, is.numeric, TRUE)
  show[num] <- lapply(show[num], function(v) round(v, 4))
  print(show, row.names = FALSE)
  if (x$order_dependent) {
    cat(sprintf("Item order changes the solution: with %s listed first the default start reaches solution %d, with %s first solution %d.\n",
                x$labels[1], x$default_first_A, x$labels[2], x$default_first_B))
  } else {
    cat("Item order does not change the default solution.\n")
  }
  nonglobal <- c(x$default_first_A, x$default_first_B)
  if (any(!is.na(nonglobal) & nonglobal > 1))
    cat("A default fit stopped at a local minimum; solution 1 has the lowest discrepancy.\n")
  invisible(x)
}

quiet_lavaan <- function(expr) {
  withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
}
