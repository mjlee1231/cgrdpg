#include <Rcpp.h>
using namespace Rcpp;

//' Fast C++ implementation of psi function
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi(s)
//' @export
// [[Rcpp::export]]
NumericVector psi_cpp(NumericVector s, double tau = 0.05) {
  int n = s.size();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double si = s[i];
    if (si < tau) {
      result[i] = tau * tau / 2.0;
    } else if (si > 1.0 - tau) {
      result[i] = si - tau / 2.0;
    } else {
      result[i] = si * si / 2.0;
    }
  }

  return result;
}

//' Fast C++ implementation of dpsi (derivative of psi)
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi'(s)
//' @export
// [[Rcpp::export]]
NumericVector dpsi_cpp(NumericVector s, double tau = 0.05) {
  int n = s.size();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double si = s[i];
    if (si < tau) {
      result[i] = 0.0;
    } else if (si > 1.0 - tau) {
      result[i] = 1.0;
    } else {
      result[i] = si;
    }
  }

  return result;
}

//' Fast C++ implementation of Psi (antiderivative of psi)
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of Psi(s)
//' @export
// [[Rcpp::export]]
NumericVector Psi_cpp(NumericVector s, double tau = 0.05) {
  int n = s.size();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double si = s[i];
    if (si < tau) {
      result[i] = tau * tau * si / 2.0;
    } else if (si > 1.0 - tau) {
      result[i] = si * si / 2.0 - tau * si / 2.0 + tau * tau / 2.0;
    } else {
      result[i] = si * si * si / 6.0;
    }
  }

  return result;
}

//' Fast C++ implementation of second derivative of psi
//'
//' @param s numeric vector of inner products
//' @param tau smoothing parameter (default 0.05)
//' @return numeric vector of psi''(s)
//' @export
// [[Rcpp::export]]
NumericVector ddpsi_cpp(NumericVector s, double tau = 0.05) {
  int n = s.size();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double si = s[i];
    if (si < tau || si > 1.0 - tau) {
      result[i] = 0.0;
    } else {
      result[i] = 1.0;
    }
  }

  return result;
}
