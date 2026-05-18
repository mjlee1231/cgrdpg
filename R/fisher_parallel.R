
#' Parallelized Fisher-scoring sweep updating X
#'
#' Parallel version of fisher_sweep_X using foreach/doParallel.
#' Updates all nodes in parallel using Jacobi-style iteration (simultaneous updates).
#'
#' @param A (n x n) symmetric adjacency matrix.
#' @param X (n x d) current embedding.
#' @param Z (p_cov x d) current Z matrix.
#' @param B (p_cov x n) covariate matrix with columns corresponding to nodes.
#' @param sign_diag (d x d) diagonal signature matrix with +1 (p) and -1 (q).
#' @param tau smoothing parameter for psi.
#' @param w_cap optional cap on psi'(s); default Inf
#' @param ls_beta line search backtracking shrink factor in (0,1); default 0.35.
#' @param ls_c Armijo line search constant in (0,1); default 1e-4 (for ascent).
#' @param ls_max maximum backtracking steps; default 30.
#' @param ncores number of cores to use; default: detectCores() - 1
#' @return Updated X.
#' @export
fisher_sweep_X_parallel <- function(
    A, X, Z, B, sign_diag,
    tau = 0.05, w_cap = Inf,
    ls_beta = 0.35, ls_c = 1e-4, ls_max = 30,
    ncores = NULL) {

  # Load parallel packages
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' required for parallel computation. Install with: install.packages('foreach')")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' required for parallel computation. Install with: install.packages('doParallel')")
  }

  n <- nrow(X); d <- ncol(X); p_cov <- nrow(Z)

  # Setup parallel backend
  if (is.null(ncores)) {
    ncores <- parallel::detectCores() - 1
  }
  ncores <- min(ncores, n, parallel::detectCores())  # Don't exceed available cores or nodes

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  # Pre-compute shared quantities (Jacobi-style: all nodes use same initial state)
  Y <- X %*% sign_diag
  Yt <- t(Y)
  ZtZ <- crossprod(Z)
  f0 <- surrogate_objective(A, X, Z, B, sign_diag, tau = tau)

  # Parallel update of all nodes
  i <- NULL  # Avoid R CMD check NOTE
  `%dopar%` <- foreach::`%dopar%`
  results <- foreach::foreach(
    i = 1:n,
    .combine = rbind,
    .packages = c("cgrdpg"),
    .export = c("dpsi", "psi", "Psi", "surrogate_objective")
  ) %dopar% {
    xi <- X[i, ]

    # Gradient pieces
    s <- as.vector(Y %*% xi)
    w <- dpsi(s, tau = tau)
    if (is.finite(w_cap)) w <- pmin(w, w_cap)
    w[i] <- 0  # No self-loops

    r <- (A[i, ] - s) * w
    s_net <- as.vector(Yt %*% r)
    bi <- B[, i]
    resid <- bi - Z %*% xi
    s_cov <- as.vector(t(Z) %*% resid)
    S <- s_net + s_cov

    # Fisher information
    G_net <- Yt %*% (Y * w)
    G_cov <- ZtZ
    G <- G_net + G_cov

    # Solve for direction
    R <- try(chol(G), silent = TRUE)
    p <- if (!inherits(R, "try-error")) {
      backsolve(R, forwardsolve(t(R), S))
    } else {
      solve(G, S)
    }

    # Check if ascent direction
    sp <- sum(S * p)
    if (!is.finite(sp) || sp <= 0) {
      p <- S
      sp <- sum(S * p)
      if (!is.finite(sp) || sp <= 0) {
        return(xi)  # Skip update
      }
    }

    # Backtracking line search
    best_eta <- 0
    eta <- 1.0

    for (bt in 1:ls_max) {
      xi_try <- xi + eta * as.vector(p)
      X_try <- X
      X_try[i, ] <- xi_try
      f_try <- surrogate_objective(A, X_try, Z, B, sign_diag, tau = tau)

      if (f_try >= f0 + ls_c * eta * sp) {
        best_eta <- eta
        break
      }
      eta <- eta * ls_beta
    }

    # Return updated position
    xi + best_eta * as.vector(p)
  }

  # Return updated X matrix
  results
}


#' Parallelized fit_grdpg_cov
#'
#' Parallel version of fit_grdpg_cov using parallelized Fisher scoring.
#'
#' @param A (n x n) symmetric adjacency matrix with entries in {0,1}.
#' @param B (p_cov x n) covariate matrix arranged as features-by-nodes.
#' @param d embedding dimension (d = p + q).
#' @param p number of positive signature directions. If \code{NULL}, estimated.
#' @param q number of negative signature directions. If \code{NULL}, taken as \code{d-p}.
#' @param tau smoothing threshold in (0, 1/2). Default 0.001.
#' @param maxit maximum number of Fisher-scoring sweeps. Default 30.
#' @param tol stopping tolerance on the maximum row change of X. Default 0.005.
#' @param ncores number of cores to use; default: detectCores() - 1
#' @return a list with \code{X}, \code{Z}, \code{Y}, \code{p}, \code{q}, \code{tau},
#' \code{converged}, \code{iters}, and \code{history}.
#' @export
fit_grdpg_cov_parallel <- function(
    A, B, d, p = NULL, q = NULL,
    tau = 0.001, maxit = 30, tol = 0.005,
    ncores = NULL
) {
  n <- nrow(A)
  p_cov <- nrow(B)

  # Input checks
  stopifnot(is.matrix(A), nrow(A) == ncol(A), is.matrix(B))
  stopifnot(ncol(B) == n)

  # Diagonal augmentation for better ASE initialization
  A_aug <- A
  diag(A_aug) <- rowSums(A) / (n - 1)

  # Initial embedding
  ase <- ase_grdpg(A_aug, d = d, p = p, q = q)
  X <- ase$X; p <- ase$p; q <- ase$q; sign_diag <- ase$sign_diag

  XtX <- crossprod(X)
  Z <- B %*% X %*% solve(XtX)

  Y <- X %*% sign_diag
  hist <- list(max_row_change = numeric(0),
               objective = numeric(0))
  converged <- FALSE

  # Initial objective
  hist$objective <- c(hist$objective,
                      surrogate_objective(A, X, Z, B, sign_diag, tau = tau))

  # PARALLEL Fisher scoring
  for (t in 1:maxit) {
    current_beta <- if (t <= 8) 0.8 else 0.4

    Xnew <- fisher_sweep_X_parallel(
      A, X, Z, B, sign_diag,
      tau = tau,
      ls_beta = current_beta,
      ls_c = 1e-4,
      ncores = ncores
    )

    # Refit Z given new X
    XtX <- crossprod(Xnew)
    Z   <- B %*% Xnew %*% solve(XtX)

    # Convergence stats
    diff <- sqrt(rowSums((Xnew - X)^2))
    mrc  <- max(diff)
    hist$max_row_change <- c(hist$max_row_change, mrc)

    # Objective
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
