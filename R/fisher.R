
#' One Fisher-scoring sweep updating X
#'
#' Implements the per-node Fisher-scoring update for the surrogate log-likelihood
#' described in the model, fixing \eqn{Y} and \eqn{Z} while updating each row \eqn{x_i}
#' with an Armijo backtracking line search (maximization).
#'
#' @param A (n x n) symmetric adjacency matrix.
#' @param X (n x d) current embedding.
#' @param Z (p_cov x d) current Z matrix.
#' @param B (p_cov x n) covariate matrix with columns corresponding to nodes.
#' @param sign_diag (d x d) diagonal signature matrix with +1 (p) and -1 (q).
#' @param tau smoothing parameter for psi.
#' @param ridge small ridge added to the information matrix for numerical stability.
#' @param w_cap optional cap on psi'(s); default Inf
#' @param ls_beta line search backtracking shrink factor in (0,1); default 0.1.
#' @param ls_c Armijo line search constant in (0,1); default 1e-2 (for ascent).
#' @param ls_max maximum backtracking steps; default 20.
#' @return Updated X.
#' @export
fisher_sweep_X <- function(
    A, X, Z, B, sign_diag,
    tau = 0.05, ridge = 0, w_cap = Inf,
    ls_beta = 0.1, ls_c = 1e-2, ls_max = 20) {

  n <- nrow(X); d <- ncol(X)

  # WORKING STATE that we keep updating as we sweep
  Xcur <- X
  Ycur <- Xcur %*% sign_diag
  Ytcur <- t(Ycur)
  ZtZ <- crossprod(Z)


  for (i in 1:n) {
    xi <- Xcur[i, ]

    # ----- gradient pieces computed at CURRENT state-----
    s <- as.vector(Ycur %*% xi)            # s_j = x_i^T y_j, n-vector
    #as.vector(xi %*% Yt)
    w <- dpsi(s, tau = tau)             # w_j = ψ'(s_j)
    if (is.finite(w_cap)) w <- pmin(w, w_cap)  # optional clipping

    r <- (A[i, ] - s) * w               # r_j = (A_ij - s_j) ψ'(s_j), n-vector
    g_net <- as.vector(Ytcur %*% r)        # Σ_j r_j y_j, d-vector
    bi <- B[, i]                        # b_i
    resid <- bi - Z %*% xi
    g_cov <- as.vector(t(Z) %*% resid)
    G <- g_net + g_cov                  # ∂l_n/∂x_i

    # ----- information and Newton/Fisher direction -----
    I_net <- Ytcur %*% (Ycur * w)             # Y^T diag(w) Y (broadcast trick)
    I_cov <- ZtZ                # I_cov <- ZtZ,
    #I <- I_net + I_cov + diag(ridge, d) # ridge adds a tiny λI for numerical stability (not in the theory; practical).
    I <- I_net + I_cov + if (ridge > 0) diag(ridge, d) else 0

    # solve I p = G for direction p (ascent) using Cholesky (stable, efficient) if I is SPD.
    # If chol fails numerically, fall back to solve().
    R <- try(chol(I), silent = TRUE)
    p <- if (!inherits(R, "try-error")) backsolve(R, forwardsolve(t(R), G)) else solve(I, G)

    # fallback: if direction is not ascent, use gradient
    gp <- sum(G * p)
    if (!is.finite(gp) || gp <= 0) {
      p  <- G
      gp <- sum(G * p)
      if (!is.finite(gp) || gp <= 0) next
    } # monitor numerically whether the Newton/Fisher direction p is truly an ascent direction
    # g^tp > 0: p is ascent (good)

    # ----- backtracking line search (Armijo for MAXIMIZATION) around CURRENT state -----
    # Condition:  ℓ(xi + η p) >= ℓ(xi) + ls_c * η * (G^T p)
    # where ℓ is the surrogate objective holding Y,Z fixed except xi.
    eta <- 1.0
    # current objective (only needs terms depending on xi, but we reuse full helper for clarity)
    f0  <- surrogate_objective(A, Xcur, Z, B, sign_diag, tau = tau)
    for (bt in 1:ls_max) {
     xi_try <- xi + eta * as.vector(p) # 1) trial update of x_i
     X_try  <- Xcur                       # 2) copy current X
     X_try[i, ] <- xi_try              #    and replace only row i
     # recompute Y for objective (only row i changes; keep it simple here)
     f_try <- surrogate_objective(A, X_try, Z, B, sign_diag, tau = tau)
     if (f_try >= f0 + ls_c * eta * gp) break # 3) Armijo (ascent) test: sufficient increase -> accept
     eta <- eta * ls_beta              # 4) otherwise shrink eta (backtrack)
    }

    # accept step and UPDATE THE WORKING STATE
    Xcur[i, ] <- xi + eta * as.vector(p)
    Ycur[i, ] <- Xcur[i, ] %*% sign_diag
    Ytcur <-  t(Ycur)
    }
  Xcur
}

#' Surrogate log-likelihood value (up to an additive constant)
#'
#' @param A (n x n) symmetric adjacency matrix.
#' @param X (n x d) embedding.
#' @param Z (p_cov x d) covariate loading matrix.
#' @param B (p_cov x n) covariate matrix (features by nodes).
#' @param sign_diag (d x d) signature matrix with +1 (p) and -1 (q).
#' @param tau smoothing parameter for \code{psi}.
#' @return Numeric scalar: the surrogate objective value.
#' @export
surrogate_objective <- function(A, X, Z, B, sign_diag, tau = 0.05) {
  # Network term: sum_{ij} [(A_ij - s_ij) ψ(s_ij) + Ψ(s_ij)], s_ij = x_i^T y_j
  Y <- X %*% sign_diag
  S <- X %*% t(Y)                              # S_{ij} = x_i^T y_j
  net <- sum((A - S) * psi(S, tau = tau) + Psi(S, tau = tau))
  # Gaussian covariate term: -1/2 ||B - Z X^T||_F^2
  cov <- -0.5 * sum((B - Z %*% t(X))^2)
  net + cov
}
