function [X_new, Z_new] = compute_cgrdpg_one_step_batch(A, B, X_init, S_estimated, tau)
% COMPUTE_CGRDPG_ONE_STEP_BATCH One-step Fisher scoring (batch update)
%
% Performs ONE Newton-Raphson update from ASE initialization.
% All vertices updated simultaneously (batch mode).
%
% Inputs:
%   A            - (n x n) adjacency matrix
%   B            - (p_cov x n) covariate matrix
%   X_init       - (n x d) initial positions (from ASE unsigned)
%   S_estimated  - (d x d) estimated signature matrix
%   tau          - clipping parameter for edge probabilities
%
% Outputs:
%   X_new  - (n x d) updated latent positions (one-step)
%   Z_new  - (p_cov x d) updated covariate loadings

n = size(A, 1);
p_cov = size(B, 1);
d = size(X_init, 2);

% Initial Z estimate: Z = B * X * (X'X)^{-1}
XtX = X_init' * X_init;
Z_init = B * (X_init / XtX);

% Batch Newton-Raphson update for all vertices
X_new = zeros(n, d);

for i = 1:n
    % Current position
    x_i = X_init(i, :)';

    % Y_i = S * x_i (for computing s_ij = x_j^T * S * x_i)
    Y_i = S_estimated * x_i;

    % Indices excluding vertex i
    idx_j = setdiff(1:n, i);

    % Edge probabilities: s_ij = x_j^T * S * x_i
    s_i = X_init(idx_j, :) * Y_i;

    % Clipped probabilities using psi function
    p_i = max(min(1 ./ (1 + exp(-s_i / tau)), 1 - 1e-10), 1e-10);

    % Residuals
    resid = A(i, idx_j)' - p_i;

    % Fisher information weights: dpsi(s, tau)
    dpsi_val = 1 ./ (tau * p_i .* (1 - p_i));

    % For gradient and Fisher info, we need Y_j = X_j * S
    Y_j = X_init(idx_j, :) * S_estimated;

    % Gradient (network part)
    grad_net = Y_j' * (resid .* dpsi_val);

    % Gradient (covariate part)
    grad_cov = Z_init' * (B(:, i) - Z_init * x_i);

    % Total gradient (scaled by n + p_cov)
    grad = (grad_net + grad_cov) / (n + p_cov);

    % Fisher information matrix (network part)
    G_net = Y_j' * (Y_j .* dpsi_val);

    % Fisher information matrix (covariate part)
    G_cov = Z_init' * Z_init;

    % Total Fisher information (scaled by n + p_cov)
    G = (G_net + G_cov) / (n + p_cov);

    % Newton-Raphson update: x_new = x_old + G^{-1} * grad
    X_new(i, :) = (x_i + (G + 1e-9 * eye(d)) \ grad)';
end

% Update Z based on new X
XtX_new = X_new' * X_new;
Z_new = B * (X_new / XtX_new);

end
