#include <Rcpp.h>
using namespace Rcpp;

//' Compute Fisher information matrix G_in for a single vertex (C++ version)
//'
//' @param i vertex index (0-based in C++)
//' @param X matrix of latent positions (n x d)
//' @param Y matrix Y = X * sign_diag (n x d)
//' @param Z covariate loading matrix (p_cov x d)
//' @param w vector of weights from dpsi (length n)
//' @return d x d Fisher information matrix
//' @export
// [[Rcpp::export]]
NumericMatrix compute_G_cpp(int i, NumericMatrix X, NumericMatrix Y,
                            NumericMatrix Z, NumericVector w) {
  int n = X.nrow();
  int d = X.ncol();
  int p_cov = Z.nrow();

  // G_net = Y^T diag(w) Y (with w[i] = 0 already set)
  NumericMatrix G_net(d, d);
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++) {
      double sum = 0.0;
      for (int j = 0; j < n; j++) {
        sum += Y(j, k) * Y(j, l) * w[j];
      }
      G_net(k, l) = sum;
    }
  }

  // G_cov = Z^T Z
  NumericMatrix G_cov(d, d);
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++) {
      double sum = 0.0;
      for (int j = 0; j < p_cov; j++) {
        sum += Z(j, k) * Z(j, l);
      }
      G_cov(k, l) = sum;
    }
  }

  // G = G_net + G_cov
  NumericMatrix G(d, d);
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++) {
      G(k, l) = G_net(k, l) + G_cov(k, l);
    }
  }

  return G;
}

//' Compute gradient (score function) for a single vertex (C++ version)
//'
//' @param i vertex index (0-based in C++)
//' @param A adjacency vector for vertex i (length n)
//' @param X matrix of latent positions (n x d)
//' @param Y matrix Y = X * sign_diag (n x d)
//' @param Z covariate loading matrix (p_cov x d)
//' @param B covariate vector for vertex i (length p_cov)
//' @param s vector of inner products (length n)
//' @param w vector of weights from dpsi (length n)
//' @return d-dimensional gradient vector
//' @export
// [[Rcpp::export]]
NumericVector compute_gradient_cpp(int i, NumericVector A, NumericMatrix X,
                                   NumericMatrix Y, NumericMatrix Z,
                                   NumericVector B, NumericVector s,
                                   NumericVector w) {
  int n = X.nrow();
  int d = X.ncol();
  int p_cov = Z.nrow();

  // Compute residuals r = (A - s) * w
  NumericVector r(n);
  for (int j = 0; j < n; j++) {
    r[j] = (A[j] - s[j]) * w[j];
  }

  // Network score: s_net = Y^T r
  NumericVector s_net(d);
  for (int k = 0; k < d; k++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum += Y(j, k) * r[j];
    }
    s_net[k] = sum;
  }

  // Get x_i
  NumericVector x_i(d);
  for (int k = 0; k < d; k++) {
    x_i[k] = X(i, k);
  }

  // Covariate residual: resid = B - Z * x_i
  NumericVector resid(p_cov);
  for (int j = 0; j < p_cov; j++) {
    double sum = 0.0;
    for (int k = 0; k < d; k++) {
      sum += Z(j, k) * x_i[k];
    }
    resid[j] = B[j] - sum;
  }

  // Covariate score: s_cov = Z^T resid
  NumericVector s_cov(d);
  for (int k = 0; k < d; k++) {
    double sum = 0.0;
    for (int j = 0; j < p_cov; j++) {
      sum += Z(j, k) * resid[j];
    }
    s_cov[k] = sum;
  }

  // Total score: S = s_net + s_cov
  NumericVector S(d);
  for (int k = 0; k < d; k++) {
    S[k] = s_net[k] + s_cov[k];
  }

  return S;
}

//' Vectorized computation of G_in matrix (optimized C++ version)
//'
//' @param i vertex index (0-based)
//' @param X matrix of latent positions (n x d)
//' @param Y matrix Y = X * sign_diag (n x d)
//' @param Z covariate loading matrix (p_cov x d)
//' @param tau smoothing parameter
//' @return d x d matrix G_in
//' @export
// [[Rcpp::export]]
NumericMatrix compute_G_in_vectorized_cpp(int i, NumericMatrix X,
                                          NumericMatrix Y, NumericMatrix Z,
                                          double tau = 0.005) {
  int n = X.nrow();
  int d = X.ncol();
  int p_cov = Z.nrow();

  // Compute s = X[i,] %*% t(Y)
  NumericVector s(n);
  for (int j = 0; j < n; j++) {
    double sum = 0.0;
    for (int k = 0; k < d; k++) {
      sum += X(i, k) * Y(j, k);
    }
    s[j] = sum;
  }

  // Compute w = dpsi(s, tau)
  NumericVector w = dpsi_cpp(s, tau);
  w[i] = 0.0;  // Exclude diagonal

  // Vectorized: G_net = t(Y * sqrt(w)) %*% (Y * sqrt(w))
  // First compute Y_weighted = Y * sqrt(w)
  NumericMatrix Y_weighted(n, d);
  for (int j = 0; j < n; j++) {
    double sqrt_w = sqrt(w[j]);
    for (int k = 0; k < d; k++) {
      Y_weighted(j, k) = Y(j, k) * sqrt_w;
    }
  }

  // G_net = crossprod(Y_weighted)
  NumericMatrix G_net(d, d);
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++) {
      double sum = 0.0;
      for (int j = 0; j < n; j++) {
        sum += Y_weighted(j, k) * Y_weighted(j, l);
      }
      G_net(k, l) = sum;
    }
  }

  // G_cov = crossprod(Z)
  NumericMatrix G_cov(d, d);
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++) {
      double sum = 0.0;
      for (int j = 0; j < p_cov; j++) {
        sum += Z(j, k) * Z(j, l);
      }
      G_cov(k, l) = sum;
    }
  }

  // G_in = (G_net + G_cov) / (n + p_cov)
  double scale = 1.0 / (n + p_cov);
  NumericMatrix G_in(d, d);
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++) {
      G_in(k, l) = (G_net(k, l) + G_cov(k, l)) * scale;
    }
  }

  return G_in;
}
