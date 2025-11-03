
#' Smoothed logx/log(1-x) functions for probabilities
#'
#' We define \code{psi1(x)} to behave like \eqn{\log x} on \eqn{[\tau, 1-\tau]}
#' and extend it quadratically to the left of \eqn{\tau} so that the function
#' and its first two derivatives are continuous. Similarly, \code{psi2(x)}
#' behaves like \eqn{\log(1-x)} on \eqn{[\tau, 1-\tau]} and is extended
#' quadratically to the right of \eqn{1-\tau}. We then set
#' \eqn{\psi(x) = \psi_1(x) - \psi_2(x)}. These are useful as a barrier-less,
#' twice-differentiable surrogate for the Bernoulli log-likelihood when dot
#' products can leave \eqn{(0,1)} during optimization.
#'
#' @param x numeric vector.
#' @param tau smoothing threshold in (0, 1/2). Default \code{0.05}.
#' @return numeric vector.
#' @examples
#' psi(0.2); dpsi(0.2); ddpsi(0.2)
#' @export
psi1 <- function(x, tau = 0.05) {
  stopifnot(tau > 0, tau < 0.5)
  a <- -1/(2*tau^2)
  b <-  2/tau
  c <-  log(tau) - 3/2
  # piecewise: quadratic if x < tau, else log(x)
  out <- ifelse(x < tau, a*x^2 + b*x + c, log(pmax(x, .Machine$double.eps)))
  return(out)
}

#' @rdname psi1
#' @export
psi2 <- function(x, tau = 0.05) {
  stopifnot(tau > 0, tau < 0.5)
  s  <- 1 - tau
  a  <- -1/(2*tau^2)
  b  <- (1 - 2*tau)/tau^2
  c  <- log(tau) + ((3*tau - 1)*(1 - tau))/(2*tau^2)  # equals log(tau) - (1/2 - 2*tau + 1.5*tau^2)/tau^2
  out <- ifelse(x > s, a*x^2 + b*x + c, log(pmax(1 - x, .Machine$double.eps)))
  return(out)
}

#' @rdname psi1
#' @export
psi <- function(x, tau = 0.05) {
  psi1(x, tau) - psi2(x, tau)
}

#' @rdname psi1
#' @export
dpsi1 <- function(x, tau = 0.05) {
  a <- -1/(2*tau^2); b <- 2/tau
  ifelse(x < tau, 2*a*x + b, 1/pmax(x, .Machine$double.eps))
}

#' @rdname psi1
#' @export
dpsi2 <- function(x, tau = 0.05) {
  s <- 1 - tau; a <- -1/(2*tau^2); b <- (1 - 2*tau)/tau^2
  ifelse(x > s, 2*a*x + b, -1/pmax(1 - x, .Machine$double.eps))
}

#' @rdname psi1
#' @export
dpsi <- function(x, tau = 0.05) {
  dpsi1(x, tau) - dpsi2(x, tau)
}

#' @rdname psi1
#' @export
ddpsi1 <- function(x, tau = 0.05) {
  a <- -1/(2*tau^2)
  ifelse(x < tau, 2*a, -1/pmax(x, .Machine$double.eps)^2)
}

#' @rdname psi1
#' @export
ddpsi2 <- function(x, tau = 0.05) {
  s <- 1 - tau; a <- -1/(2*tau^2)
  ifelse(x > s, 2*a, -1/pmax(1 - x, .Machine$double.eps)^2)
}

#' @rdname psi1
#' @export
ddpsi <- function(x, tau = 0.05) {
  ddpsi1(x, tau) - ddpsi2(x, tau)
}

#' Primitive (integral) of psi1
#' @keywords internal
Psi1 <- function(x, tau = 0.05) {
  stopifnot(tau > 0, tau < 0.5)
  a <- -1/(2*tau^2); b <- 2/tau; c <- log(tau) - 3/2
  # integral of quadratic from 0 to t: a/3 t^3 + b/2 t^2 + c t
  Q <- function(t) a/3*t^3 + b/2*t^2 + c*t
  out <- numeric(length(x))
  # constant so that piecewise is continuous at tau
  Ctau <- Q(tau)
  idx_q <- (x < tau)
  out[idx_q] <- Q(x[idx_q])
  # integral of log from tau to x: x log x - x - (tau log tau - tau)
  out[!idx_q] <- Ctau + (x[!idx_q]*log(pmax(x[!idx_q], .Machine$double.eps)) - x[!idx_q]) -
    (tau*log(tau) - tau)
  out
}

#' Primitive (integral) of psi2
#' @keywords internal
Psi2 <- function(x, tau = 0.05) {
  stopifnot(tau > 0, tau < 0.5)
  s <- 1 - tau
  a <- -1/(2*tau^2); b <- (1 - 2*tau)/tau^2
  c <- log(tau) + ((3*tau - 1)*(1 - tau))/(2*tau^2)
  # quadratic primitive
  Q <- function(t) a/3*t^3 + b/2*t^2 + c*t
  out <- numeric(length(x))
  idx_log <- (x <= s)
  # ∫_0^x log(1-t) dt = -(1-x) log(1-x) - x
  out[idx_log] <- -(1 - x[idx_log]) * log(pmax(1 - x[idx_log], .Machine$double.eps)) - x[idx_log]
  if (any(!idx_log)) {
    # constant: integral of log(1-t) from 0 to s
    Clog <- -(1 - s) * log(tau) - s   # since 1 - s = tau
    out[!idx_log] <- Clog + (Q(x[!idx_log]) - Q(s))
  }
  out
}

#' Primitive (integral) of psi = psi1 - psi2
#' @rdname psi1
#' @export
Psi <- function(x, tau = 0.05) {
  Psi1(x, tau) - Psi2(x, tau)
}

