function [f, grad, H] = surrogate_objective_2d(x, A, B, n, p_cov, d, tau)
% Surrogate objective for 2D RDPG with covariates (Z0 = 0, B = noise)
% Uses pseudo-likelihood with psi surrogate for logistic link
%
% Inputs:
%   x      - stacked vector [X(:); Y(:); Z(:)]
%   A      - (n x n) adjacency matrix
%   B      - (p_cov x n) covariate matrix (pure noise)
%   n      - number of vertices
%   p_cov  - number of covariates
%   d      - embedding dimension (2)
%   tau    - surrogate parameter
%
% Outputs:
%   f    - objective value
%   grad - gradient vector
%   H    - Hessian matrix

% Extract variables
X = reshape(x(1:n*d), n, d);
Y = reshape(x(n*d+1:2*n*d), n, d);
Z = reshape(x(2*n*d+1:end), p_cov, d);

% Network term: psi surrogate for A
S_net = X * Y';  % (n x n)
[psi_vals, dpsi_vals, ddpsi_vals] = psi_all(S_net(:), tau);
Psi_net = reshape(psi_vals, n, n);
dPsi_net = reshape(dpsi_vals, n, n);
ddPsi_net = reshape(ddpsi_vals, n, n);

% Network pseudo-likelihood
f_net = -sum(A(:) .* S_net(:)) + sum(Psi_net(:));

% Regression term: ||B - Z*X'||_F^2
Resid = B - Z * X';
f_reg = 0.5 * sum(Resid(:).^2);

% Total objective
f = f_net + f_reg;

if nargout > 1
    % Gradient computation
    grad_X = zeros(n, d);
    grad_Y = zeros(n, d);
    grad_Z = zeros(p_cov, d);

    % Network gradient
    % ∂f_net/∂X = Y * (dPsi - A)'
    % ∂f_net/∂Y = X' * (dPsi - A)
    Diff_net = dPsi_net - A;
    grad_X = grad_X + Diff_net * Y;
    grad_Y = grad_Y + Diff_net' * X;

    % Regression gradient
    % ∂f_reg/∂X = -Z' * (B - Z*X')
    % ∂f_reg/∂Z = -(B - Z*X') * X
    grad_X = grad_X - Z' * Resid;
    grad_Z = grad_Z - Resid * X;

    grad = [grad_X(:); grad_Y(:); grad_Z(:)];
end

if nargout > 2
    % Hessian computation (block diagonal approximation)
    % Network Hessian blocks
    W_net = ddPsi_net;  % (n x n)

    % H_XX ≈ Y' * diag(sum(W, 2)) * Y
    % H_YY ≈ X' * diag(sum(W, 1)) * X
    W_row_sums = sum(W_net, 2);
    W_col_sums = sum(W_net, 1)';

    H_XX = kron(eye(d), diag(W_row_sums)) + kron(Z' * Z, eye(n));
    H_YY = kron(eye(d), diag(W_col_sums));
    H_ZZ = kron(X * X', eye(p_cov));

    % Block diagonal Hessian
    H = blkdiag(H_XX, H_YY, H_ZZ);
end

end

function [psi, dpsi, ddpsi] = psi_all(s, tau)
% Compute psi, first derivative, and second derivative
% psi(s; tau) = log(1 + exp(s)) + tau * s^2

psi = log(1 + exp(s)) + tau * s.^2;

if nargout > 1
    sig = 1 ./ (1 + exp(-s));
    dpsi = sig + 2 * tau * s;
end

if nargout > 2
    ddpsi = sig .* (1 - sig) + 2 * tau;
end
end
