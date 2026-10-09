#' Population correlation matrix of a two-trait design
#'
#' Items in block A load only on trait A and items in block B only on trait
#' B; the traits correlate `phi`, and every item has unit variance.
#'
#' @param phi Trait correlation, strictly between -1 and 1.
#' @param lambda_A Standardized loadings of block A: one value shared by all
#'   `m_A` items, or one value per item.
#' @param lambda_B Loadings of block B, in the same form. Defaults to
#'   `lambda_A`.
#' @param m_A,m_B Numbers of items, used when the loadings are given as single
#'   values. Ignored for a block whose loadings are given per item.
#' @return The `(m_A + m_B)` by `(m_A + m_B)` correlation matrix, with items
#'   named `A1, A2, ..., B1, B2, ...`. The design is attached as attribute
#'   `"two_trait"` and is used by [pop_fit()] to compute the share of each
#'   trait.
#' @examples
#' two_trait_cor(phi = 0.5, lambda_A = 0.7, m_A = 4)
#' two_trait_cor(phi = 0.3, lambda_A = c(0.6, 0.7, 0.8), lambda_B = 0.75, m_B = 5)
#' @export
two_trait_cor <- function(phi, lambda_A, lambda_B = lambda_A, m_A = 4, m_B = m_A) {
  check_phi(phi)
  if (length(phi) != 1) stop("`phi` must be a single value.", call. = FALSE)
  check_loading(lambda_A, "lambda_A")
  check_loading(lambda_B, "lambda_B")
  if (length(lambda_A) > 1) m_A <- length(lambda_A)
  if (length(lambda_B) > 1) m_B <- length(lambda_B)
  check_m(m_A, "m_A"); check_m(m_B, "m_B")
  lA <- rep_len(lambda_A, m_A); lB <- rep_len(lambda_B, m_B)
  p <- m_A + m_B
  Lam <- matrix(0, p, 2)
  Lam[seq_len(m_A), 1] <- lA
  Lam[m_A + seq_len(m_B), 2] <- lB
  R <- Lam %*% matrix(c(1, phi, phi, 1), 2) %*% t(Lam)
  diag(R) <- 1
  nm <- c(paste0("A", seq_len(m_A)), paste0("B", seq_len(m_B)))
  dimnames(R) <- list(nm, nm)
  attr(R, "two_trait") <- list(phi = phi, lambda_A = lA, lambda_B = lB,
                               blocks = list(A = seq_len(m_A), B = m_A + seq_len(m_B)))
  R
}
