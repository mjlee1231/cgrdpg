
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
#' @param w_cap optional cap on psi'(s); default Inf
#' @param ls_beta line search backtracking shrink factor in (0,1); default 0.35.
#' @param ls_c Armijo line search constant in (0,1); default 1e-4 (for ascent).
#' @param ls_max maximum backtracking steps; default 30.
#' @return Updated X.
#' @export
fisher_sweep_X <- function(
    A, X, Z, B, sign_diag,
    tau = 0.05, w_cap = Inf,
    ls_beta = 0.35, ls_c = 1e-4, ls_max = 30) {

  n <- nrow(X); d <- ncol(X); p_cov <- nrow(Z)

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

    # IMPORTANT: For graphs without self-loops (diag(A)=0), exclude j=i term
    w[i] <- 0

    r <- (A[i, ] - s) * w               # r_j = (A_ij - s_j) ψ'(s_j), n-vector (j≠i)
    s_net <- as.vector(Ytcur %*% r)        # Σ_{j≠i} r_j y_j, d-vector
    bi <- B[, i]                        # b_i
    resid <- bi - Z %*% xi
    s_cov <- as.vector(t(Z) %*% resid)
    S <- s_net + s_cov  # score function

    # ----- information and Newton/Fisher direction -----
    # Fisher information for node i: sum over all edges incident to i
    # For graphs without self-loops, sum over j≠i only (w[i] already set to 0)
    # MODIFICATION: Truncate s to [tau, 1-tau] in Fisher info denominator only
    s_clipped <- pmax(pmin(s, 1 - tau), tau)
    w_fisher <- dpsi(s_clipped, tau = tau)
    if (is.finite(w_cap)) w_fisher <- pmin(w_fisher, w_cap)
    w_fisher[i] <- 0
    G_net <- Ytcur %*% (Ycur * w_fisher)      # Y^T diag(w_fisher) Y, with w[i]=0
    G_cov <- ZtZ                              # G_cov <- ZtZ
    G <- G_net + G_cov        # Total Fisher information 

    # solve G p = S for direction p (ascent) using Cholesky (stable, efficient) if G is SPD.
    # If chol fails numerically, fall back to solve().
    R <- try(chol(G), silent = TRUE)
    p <- if (!inherits(R, "try-error")) backsolve(R, forwardsolve(t(R), S)) else solve(G, S)

    # fallback: if direction is not ascent, use gradient
    # sp <- sum(S * p)
    # if (!is.finite(sp) || sp <= 0) {
    #   p  <- S
    #   sp <- sum(S * p)
    #   if (!is.finite(sp) || sp <= 0) next
    # } # monitor numerically whether the Newton/Fisher direction p is truly an ascent direction
    # # s^tp > 0: p is ascent (good)
    
    # ----- FALLBACK MONITORING -----
    sp <- sum(S * p)
    if (!is.finite(sp) || sp <= 0) {
      message(sprintf("Warning: Node %d - Fisher direction is not ascent (sp = %.2e). Falling back to Gradient.", i, sp))
      p  <- S
      sp <- sum(S * p)
      if (!is.finite(sp) || sp <= 0) {
        message(sprintf("Critical: Node %d - Gradient itself is invalid. Skipping update.", i))
        next
      }
    }

    # ----- backtracking line search (Armijo for MAXIMIZATION) around CURRENT state -----
    # Condition:  l(xi + η p) >= l(xi) + ls_c * η * (S^T p)
    # where l is the surrogate objective holding Y,Z fixed except xi.
    best_eta <- 0 # basically not changing
    eta <- 1.0
    f0 <- surrogate_objective(A, Xcur, Z, B, sign_diag, tau = tau)
    
    for (bt in 1:ls_max) {
      xi_try <- xi + eta * as.vector(p)
      X_try <- Xcur
      X_try[i, ] <- xi_try
      f_try <- surrogate_objective(A, X_try, Z, B, sign_diag, tau = tau)
      
      if (f_try >= f0 + ls_c * eta * sp) {
        best_eta <- eta # save successful eta
        break 
      }
      eta <- eta * ls_beta
    }

    # accept step and UPDATE THE WORKING STATE
    Xcur[i, ] <- xi + best_eta * as.vector(p)
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
  # Network term: sum over (i,j) pairs with i≠j: Σ_i Σ_{j≠i} [(A_ij - s_ij) ψ(s_ij) + Ψ(s_ij)]
  # For graphs without self-loops (diag(A)=0), exclude diagonal terms
  Y <- X %*% sign_diag
  S <- X %*% t(Y)                              # S_{ij} = x_i^T y_j
  diag(S) <- 0                                 # Exclude diagonal (no self-loops)
  net <- sum((A - S) * psi(S, tau = tau) + Psi(S, tau = tau))
  # Gaussian covariate term: -1/2 ||B - Z X^T||_F^2
  cov <- -0.5 * sum((B - Z %*% t(X))^2)
  net + cov
}
