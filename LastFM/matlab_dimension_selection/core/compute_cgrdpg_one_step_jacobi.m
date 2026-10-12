function [X_new, Z_new] = compute_cgrdpg_one_step_jacobi(A, B, X_init, S_estimated, tau, step_size)
% COMPUTE_CGRDPG_ONE_STEP_JACOBI One-step Fisher scoring (Jacobi update)
%
% Performs ONE Newton-Raphson update from ASE initialization.
% All vertices updated in parallel (Jacobi/vectorized mode).
%
% Inputs:
%   A            - (n x n) adjacency matrix
%   B            - (p_cov x n) covariate matrix
%   X_init       - (n x d) initial positions (from ASE unsigned)
%   S_estimated  - (d x d) estimated signature matrix
%   tau          - clipping parameter for edge probabilities
%   step_size    - step size for Jacobi update (default: 1.0 = full Newton)
%
% Outputs:
%   X_new  - (n x d) updated latent positions (one-step)
%   Z_new  - (p_cov x d) updated covariate loadings

if nargin < 6
    step_size = 1.0;  % Full Newton step
end

n = size(A, 1);
p_cov = size(B, 1);
d = size(X_init, 2);

% Initial Z estimate
XtX = X_init' * X_init;
Z_init = B * (X_init / XtX);

% Vectorized Jacobi update: compute all updates simultaneously
% then apply all at once

% Initialize arrays for vectorized computation
X_updates = zeros(n, d);

for i = 1:n
    % Current position
    x_i = X_init(i, :)';

    % Y_i = S * x_i (for computing s_ij = x_j^T * S * x_i)
    Y_i = S_estimated * x_i;

    % Indices excluding vertex i
    idx_j = setdiff(1:n, i);

    % Linear predictor: s_ij = x_j^T * S * x_i
    s_i = X_init(idx_j, :) * Y_i;

    % Residuals (using psi link, s_i is treated as pseudo-probability)
    resid = A(i, idx_j)' - s_i;

    % Fisher weights: dpsi(s, tau) using psi link
    dpsi_val = dpsi(s_i, tau);

    % For gradient and Fisher info, we need Y_j = X_j * S
    Y_j = X_init(idx_j, :) * S_estimated;

    % Gradient
    grad_net = Y_j' * (resid .* dpsi_val);
    grad_cov = Z_init' * (B(:, i) - Z_init * x_i);
    grad = (grad_net + grad_cov) / (n + p_cov);

    % Fisher information
    G_net = Y_j' * (Y_j .* dpsi_val);
    G_cov = Z_init' * Z_init;
    G = (G_net + G_cov) / (n + p_cov);

    % Compute update direction (but don't apply yet)
    X_updates(i, :) = ((G + 1e-9 * eye(d)) \ grad)';
end

% Apply all updates simultaneously (Jacobi style)
X_new = X_init + step_size * X_updates;

% Update Z based on new X
XtX_new = X_new' * X_new;
Z_new = B * (X_new / XtX_new);

end
