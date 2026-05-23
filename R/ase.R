
#' Adjacency Spectral Embedding (GRDPG)
#'
#' Compute the \eqn{d}-dimensional adjacency spectral embedding of a symmetric adjacency matrix.
#' We order eigenpairs by descending absolute eigenvalues (appropriate for GRDPG) and form
#' \eqn{X = U |\Lambda|^{1/2}} (unsigned) and \eqn{X_signed = U |\Lambda|^{1/2} sign(\Lambda)} (signed).
#' If the signature \code{(p,q)} is not supplied, it is estimated by counting the number of
#' positive eigenvalues among the selected \eqn{d}.
#'
#' @param A numeric symmetric matrix (n x n) adjacency.
#' @param d embedding dimension.
#' @param p number of positive eigenvalues. If \code{NULL}, estimated.
#' @param q number of negative eigenvalues. If \code{NULL}, taken as \code{d-p}.
#' @return a list with:
#'   \itemize{
#'     \item \code{X}: unsigned embedding \eqn{U |\Lambda|^{1/2}} (for RDPG with S=I)
#'     \item \code{X_signed}: signed embedding \eqn{U |\Lambda|^{1/2} sign(\Lambda)} (for GRDPG with indefinite S)
#'     \item \code{eigvals}: eigenvalues
#'     \item \code{U}: eigenvectors
#'     \item \code{p}, \code{q}: signature parameters
#'     \item \code{sign_diag}: signature matrix
#'   }
#' @export
ase_grdpg <- function(A, d, p = NULL, q = NULL) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  # Use RSpectra::eigs_sym for faster computation (only computes top d eigenpairs)
  eigenA <- RSpectra::eigs_sym(A, k = d, which = "LM")  # "LM" = largest magnitude (absolute value)
  vals <- eigenA$values
  U <- eigenA$vectors

  # Unsigned version: X = U * |Lambda|^{1/2}
  # Use for RDPG (S = I, all eigenvalues positive)
  X <- U %*% diag(sqrt(abs(vals)), nrow = d, ncol = d)

  # Signed version: X_signed = U * |Lambda|^{1/2} * sign(Lambda)
  # Use for GRDPG (indefinite S with negative eigenvalues)
  # This incorporates the signature into the embedding itself
  X_signed <- U %*% diag(sqrt(abs(vals)) * sign(vals), nrow = d, ncol = d)

  if (is.null(p)) {
    p <- sum(vals > 0)
  }
  if (is.null(q)) q <- d - p
  stopifnot(p + q == d)
  sign_diag <- diag(c(rep(1, p), rep(-1, q)))

  list(X = X, X_signed = X_signed, eigvals = vals, U = U, p = p, q = q, sign_diag = sign_diag)
}
