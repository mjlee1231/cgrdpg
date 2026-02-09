#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Fast C++ implementation of psi1 function
//'
//' Smoothed log(x) function: quadratic for x < tau, log(x) otherwise
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi1(x)
//' @export
// [[Rcpp::export]]
NumericVector psi1_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double a = -1.0 / (2.0 * tau * tau);
  double b = 2.0 / tau;
  double c = log(tau) - 1.5;

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi < tau) {
      result[i] = a * xi * xi + b * xi + c;
    } else {
      result[i] = log(std::max(xi, 2.220446e-16));  // .Machine$double.eps in R
    }
  }

  return result;
}

//' Fast C++ implementation of psi2 function
//'
//' Smoothed log(1-x) function: log(1-x) for x <= 1-tau, quadratic otherwise
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi2(x)
//' @export
// [[Rcpp::export]]
NumericVector psi2_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double s = 1.0 - tau;
  double a = -1.0 / (2.0 * tau * tau);
  double b = (1.0 - 2.0 * tau) / (tau * tau);
  double c = log(tau) + ((3.0 * tau - 1.0) * (1.0 - tau)) / (2.0 * tau * tau);

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi > s) {
      result[i] = a * xi * xi + b * xi + c;
    } else {
      result[i] = log(std::max(1.0 - xi, 2.220446e-16));
    }
  }

  return result;
}

//' Fast C++ implementation of psi function
//'
//' psi(x) = psi1(x) - psi2(x)
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi(s)
//' @export
// [[Rcpp::export]]
NumericVector psi_cpp(NumericVector s, double tau = 0.05) {
  return psi1_cpp(s, tau) - psi2_cpp(s, tau);
}

//' Fast C++ implementation of dpsi1 (derivative of psi1)
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi1'(x)
//' @export
// [[Rcpp::export]]
NumericVector dpsi1_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double a = -1.0 / (2.0 * tau * tau);
  double b = 2.0 / tau;

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi < tau) {
      result[i] = 2.0 * a * xi + b;
    } else {
      result[i] = 1.0 / std::max(xi, 2.220446e-16);
    }
  }

  return result;
}

//' Fast C++ implementation of dpsi2 (derivative of psi2)
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi2'(x)
//' @export
// [[Rcpp::export]]
NumericVector dpsi2_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double s = 1.0 - tau;
  double a = -1.0 / (2.0 * tau * tau);
  double b = (1.0 - 2.0 * tau) / (tau * tau);

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi > s) {
      result[i] = 2.0 * a * xi + b;
    } else {
      result[i] = -1.0 / std::max(1.0 - xi, 2.220446e-16);
    }
  }

  return result;
}

//' Fast C++ implementation of dpsi (derivative of psi)
//'
//' dpsi(x) = dpsi1(x) - dpsi2(x)
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi'(s)
//' @export
// [[Rcpp::export]]
NumericVector dpsi_cpp(NumericVector s, double tau = 0.05) {
  return dpsi1_cpp(s, tau) - dpsi2_cpp(s, tau);
}

//' Fast C++ implementation of ddpsi1 (second derivative of psi1)
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi1''(x)
//' @export
// [[Rcpp::export]]
NumericVector ddpsi1_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double a = -1.0 / (2.0 * tau * tau);

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi < tau) {
      result[i] = 2.0 * a;
    } else {
      double xi_safe = std::max(xi, 2.220446e-16);
      result[i] = -1.0 / (xi_safe * xi_safe);
    }
  }

  return result;
}

//' Fast C++ implementation of ddpsi2 (second derivative of psi2)
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi2''(x)
//' @export
// [[Rcpp::export]]
NumericVector ddpsi2_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double s = 1.0 - tau;
  double a = -1.0 / (2.0 * tau * tau);

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi > s) {
      result[i] = 2.0 * a;
    } else {
      double one_minus_xi = std::max(1.0 - xi, 2.220446e-16);
      result[i] = -1.0 / (one_minus_xi * one_minus_xi);
    }
  }

  return result;
}

//' Fast C++ implementation of ddpsi (second derivative of psi)
//'
//' ddpsi(x) = ddpsi1(x) - ddpsi2(x)
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi''(s)
//' @export
// [[Rcpp::export]]
NumericVector ddpsi_cpp(NumericVector s, double tau = 0.05) {
  return ddpsi1_cpp(s, tau) - ddpsi2_cpp(s, tau);
}

//' Fast C++ implementation of Psi1 (integral of psi1)
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of Psi1(x)
//' @keywords internal
// [[Rcpp::export]]
NumericVector Psi1_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double a = -1.0 / (2.0 * tau * tau);
  double b = 2.0 / tau;
  double c = log(tau) - 1.5;

  // Q(t) = a/3 * t^3 + b/2 * t^2 + c * t
  auto Q = [a, b, c](double t) {
    return a / 3.0 * t * t * t + b / 2.0 * t * t + c * t;
  };

  double C_tau = Q(tau);

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi < tau) {
      result[i] = Q(xi);
    } else {
      double xi_safe = std::max(xi, 2.220446e-16);
      // integral of log from tau to x: x*log(x) - x - (tau*log(tau) - tau)
      result[i] = C_tau + (xi_safe * log(xi_safe) - xi_safe) - (tau * log(tau) - tau);
    }
  }

  return result;
}

//' Fast C++ implementation of Psi2 (integral of psi2)
//'
//' @param x numeric vector
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of Psi2(x)
//' @keywords internal
// [[Rcpp::export]]
NumericVector Psi2_cpp(NumericVector x, double tau = 0.05) {
  int n = x.size();
  NumericVector result(n);

  double s = 1.0 - tau;
  double a = -1.0 / (2.0 * tau * tau);
  double b = (1.0 - 2.0 * tau) / (tau * tau);
  double c = log(tau) + ((3.0 * tau - 1.0) * (1.0 - tau)) / (2.0 * tau * tau);

  auto Q = [a, b, c](double t) {
    return a / 3.0 * t * t * t + b / 2.0 * t * t + c * t;
  };

  // integral of log(1-t) from 0 to s
  double C_log = -(1.0 - s) * log(tau) - s;  // since 1 - s = tau

  for (int i = 0; i < n; i++) {
    double xi = x[i];
    if (xi <= s) {
      // ∫_0^x log(1-t) dt = -(1-x) * log(1-x) - x
      double one_minus_xi = std::max(1.0 - xi, 2.220446e-16);
      result[i] = -(1.0 - xi) * log(one_minus_xi) - xi;
    } else {
      result[i] = C_log + (Q(xi) - Q(s));
    }
  }

  return result;
}

//' Fast C++ implementation of Psi (integral of psi)
//'
//' Psi(x) = Psi1(x) - Psi2(x)
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of Psi(s)
//' @export
// [[Rcpp::export]]
NumericVector Psi_cpp(NumericVector s, double tau = 0.05) {
  return Psi1_cpp(s, tau) - Psi2_cpp(s, tau);
}
