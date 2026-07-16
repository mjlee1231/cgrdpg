function G = compute_fisher_info_cgrdpg(i, X_mat, Y_mat, Z_mat, tau)
% COMPUTE_FISHER_INFO_CGRDPG Compute Fisher information matrix for vertex i
%
% Computes the Fisher information matrix used for inference in cgrdpg model.
% This is the precision matrix (inverse covariance) for the asymptotic
% normal distribution of the latent position estimate.
%
% Inputs:
%   i     - vertex index
%   X_mat - (n x d) latent positions
%   Y_mat - (n x d) signed latent positions (Y = X * S)
%   Z_mat - (p_cov x d) covariate coefficient matrix
%   tau   - smoothing parameter for psi functions
%
% Outputs:
%   G - (d x d) Fisher information matrix, normalized by (n + p_cov)
%
% Formula:
%   G = (Y^T * diag(dpsi(s_i)) * Y + Z^T * Z) / (n + p_cov)
% where s_i = X[i,:] * Y^T (edge probabilities from vertex i)

n = size(X_mat, 1);
p_cov = size(Z_mat, 1);

% Compute edge probabilities from vertex i to all others
s_i = X_mat(i, :) * Y_mat';  % (1 x n) vector

% Compute weights dpsi(s_i)
w = dpsi(s_i(:), tau);  % (n x 1) vector
w(i) = 0;  % Exclude self-loop

% Network Fisher information: Y^T * diag(sqrt(w)) * Y
% Note: Using sqrt(w) to match R implementation
G_net = Y_mat' * (Y_mat .* sqrt(w));

% Covariate Fisher information: Z^T * Z
G_cov = Z_mat' * Z_mat;

% Total Fisher information (normalized)
G = (G_net + G_cov) / (n + p_cov);

end
