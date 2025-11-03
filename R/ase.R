
#' Adjacency Spectral Embedding (GRDPG)
#'
#' Compute the \eqn{d}-dimensional adjacency spectral embedding of a symmetric adjacency matrix.
#' We order eigenpairs by descending absolute eigenvalues (appropriate for GRDPG) and form
#' \eqn{X = U |\Lambda|^{1/2}}. If the signature \code{(p,q)} is not supplied, it is estimated
#' by counting the number of positive eigenvalues among the selected \eqn{d}.
#'
#' @param A numeric symmetric matrix (n x n) adjacency.
#' @param d embedding dimension.
#' @param p number of positive eigenvalues. If \code{NULL}, estimated.
#' @param q number of negative eigenvalues. If \code{NULL}, taken as \code{d-p}.
#' @return a list with \code{X}, \code{eigvals}, \code{U}, \code{p}, \code{q}, \code{sign_diag}.
#' @export
ase_grdpg <- function(A, d, p = NULL, q = NULL) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  ev <- eigen(A, symmetric = TRUE)
  ord <- order(abs(ev$values), decreasing = TRUE)
  vals <- ev$values[ord][1:d]
  U <- ev$vectors[, ord, drop = FALSE][, 1:d, drop = FALSE]
  X <- U %*% diag(sqrt(abs(vals)), nrow = d, ncol = d)
  if (is.null(p)) {
    p <- sum(vals > 0)
  }
  if (is.null(q)) q <- d - p
  stopifnot(p + q == d)
  sign_diag <- diag(c(rep(1, p), rep(-1, q)))
  list(X = X, eigvals = vals, U = U, p = p, q = q, sign_diag = sign_diag)
}
