function [X_aligned, Q] = procrustes_align(X_est, X_target)
% PROCRUSTES_ALIGN Align estimated latent positions to target positions
%
% Finds orthogonal matrix Q that minimizes ||X_est * Q - X_target||_F
% using the Procrustes solution.
%
% Inputs:
%   X_est    - (n x d) estimated latent positions
%   X_target - (n x d) target (true) latent positions
%
% Outputs:
%   X_aligned - (n x d) aligned positions = X_est * Q
%   Q         - (d x d) orthogonal rotation matrix
%
% Algorithm:
%   1. Compute SVD of X_est' * X_target = U * Sigma * V'
%   2. Optimal rotation: Q = U * V'
%   3. Aligned positions: X_aligned = X_est * Q

% Compute cross-product matrix
M = X_est' * X_target;

% SVD
[U, ~, V] = svd(M);

% Optimal orthogonal transformation
Q = U * V';

% Apply transformation
X_aligned = X_est * Q;

end
