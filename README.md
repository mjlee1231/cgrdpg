
# cgrdpg

A minimal R package that implements a **generalized random dot product graph (GRDPG)** estimator with **contextual covariates**:
\[
A_{ij} \sim \mathrm{Bernoulli}(x_i^\top I_{p,q} x_j), \qquad B = Z X^\top + W,\; W_{li} \stackrel{iid}{\sim} \mathcal{N}(0,1).
\]

The algorithm uses:
1. **Adjacency Spectral Embedding (ASE)** with absolute eigenvalue ordering (GRDPG).
2. **Smoothed** surrogate log-likelihood via \(\psi\)-functions that match \(\log x\) and \(\log(1-x)\) on \([\tau, 1-\tau]\).
3. A **Fisher scoring** sweep to update \(X\) with \(Y = X I_{p,q}\) and current \(Z\) fixed.
4. **Closed-form refit** of \(Z\) each iteration: \( \hat Z = B X (X^\top X + \lambda I)^{-1} \).



## Install (devtools)

```r
install.packages("devtools")
devtools::document("cgrdpg")   # generate NAMESPACE + man/
devtools::install("cgrdpg")
```

## Example

```r
library(cgrdpg)
set.seed(1)
n <- 60; d <- 3; p_sig <- 2; q_sig <- 1; p_cov <- 4
X0 <- matrix(rnorm(n*d, sd = 0.7), n, d)
S  <- diag(c(rep(1, p_sig), rep(-1, q_sig)))
P <- X0 %*% S %*% t(X0); P <- pmin(pmax(P, 1e-3), 1-1e-3)
A <- matrix(rbinom(n*n, 1, as.vector(P)), n, n); A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0
Z0 <- matrix(rnorm(p_cov*d), p_cov, d)
B <- Z0 %*% t(X0) + matrix(rnorm(p_cov*n), p_cov, n)

fit <- fit_grdpg_cov(A, B, d = d, p = p_sig, q = q_sig, maxit = 3)
str(fit)
```

### Key knobs

- `tau`: smoothing threshold (default 1e-2). Smaller = closer to raw log-likelihood but less stable.
- `alpha`: weight on the covariate term in the Fisher update.
---
