
#' Fit GRDPG with contextual covariates via Fisher scoring
#'
#' @param A (n x n) symmetric adjacency matrix with entries in {0,1}.
#' @param B (p_cov x n) covariate matrix arranged as features-by-nodes.
#' @param d embedding dimension (d = p + q).
#' @param p number of positive signature directions. If \code{NULL}, estimated.
#' @param q number of negative signature directions. If \code{NULL}, taken as \code{d-p}.
#' @details Dimension consistency (p + q == d) is validated internally by ase_grdpg().
#' @param tau smoothing threshold in (0, 1/2). Default 0.05.
#' @param maxit maximum number of Fisher-scoring sweeps. Default 5.
#' @param tol stopping tolerance on the maximum row change of X. Default 1e-2.
#' @param ridge small ridge inside linear solves. Default 0.
#' @return a list with \code{X}, \code{Z}, \code{Y}, \code{p}, \code{q}, \code{tau},
#' \code{converged}, \code{iters}, and \code{history}.
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 60; d <- 3; p_sig <- 2; q_sig <- 1; p_cov <- 4
#' X0 <- matrix(rnorm(n*d, sd = 0.7), n, d)
#' S <- diag(c(rep(1, p_sig), rep(-1, q_sig)))
#' P <- X0 %*% S %*% t(X0)
#' A <- matrix(rbinom(n*n, 1, as.vector(P)), n, n); A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0
#' Z0 <- matrix(rnorm(p_cov*d), p_cov, d)
#' B <- Z0 %*% t(X0) + matrix(rnorm(p_cov*n), p_cov, n)
#' fit <- fit_grdpg_cov(A, B, d = d, p = p_sig, q = q_sig, maxit = 3)
#' str(fit)
#' }
#' @export
fit_grdpg_cov <- function(
    A, B, d, p = NULL, q = NULL,
    tau = 0.05, maxit = 5, tol = 1e-2, ridge = 0
) {
  # Define n first
  n <- nrow(A)
  p_cov <- nrow(B)
  
  # Input checks
  stopifnot(is.matrix(A), nrow(A) == ncol(A), is.matrix(B))
  stopifnot(ncol(B) == n)

  ase <- ase_grdpg(A, d = d, p = p, q = q)
  X <- ase$X; p <- ase$p; q <- ase$q; sign_diag <- ase$sign_diag

  XtX <- crossprod(X)
  Z <- B %*% X %*% solve(XtX)

  Y <- X %*% sign_diag
  hist <- list(max_row_change = numeric(0),
               objective = numeric(0))
  converged <- FALSE

  # initial objective
  hist$objective <- c(hist$objective,
                      surrogate_objective(A, X, Z, B, sign_diag, tau = tau))

  for (t in 1:maxit) {
    Xnew <- fisher_sweep_X(A, X, Z, B, sign_diag, tau = tau, ridge = ridge)

    # refit Z given new X
    XtX <- crossprod(Xnew)
    Z   <- B %*% Xnew %*% solve(XtX)

    # convergence stats
    diff <- sqrt(rowSums((Xnew - X)^2))
    mrc  <- max(diff)
    hist$max_row_change <- c(hist$max_row_change, mrc)

    # objective
    obj <- surrogate_objective(A, Xnew, Z, B, sign_diag, tau = tau)
    hist$objective <- c(hist$objective, obj)

    X <- Xnew
    Y <- X %*% sign_diag
    if (mrc < tol) { converged <- TRUE; break }
  }

  list(
    X = X, Z = Z, Y = Y,
    p = p, q = q, tau = tau,
    converged = converged,
    iters = t,
    history = hist
  )
}
