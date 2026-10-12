function Prec = compute_prec_ase(i, X_mat, S, eps_clip)
% COMPUTE_PREC_ASE Compute precision matrix for ASE inference
%
% Matches R's compute_prec_ase function
% Formula: Prec = S %*% Delta %*% M^{-1} %*% Delta %*% S
% where Delta = X^T X and M = X_{-i}^T diag(p*(1-p)) X_{-i}
%
% Inputs:
%   i        - vertex index
%   X_mat    - (n x d) latent positions
%   S        - (d x d) signature matrix
%   eps_clip - clipping value for probabilities
%
% Outputs:
%   Prec - (d x d) precision matrix

n = size(X_mat, 1);
d = size(X_mat, 2);

% All vertices except i
idx_j = setdiff(1:n, i);

% Edge probabilities: p_ij = x_i^T S x_j
p_vals = X_mat(idx_j, :) * (S * X_mat(i, :)');
p_vals = max(min(p_vals, 1 - eps_clip), eps_clip);

% Delta = X^T X (second moment matrix)
Delta = X_mat' * X_mat;

% M = X_{-i}^T diag(p*(1-p)) X_{-i}
weights = p_vals .* (1 - p_vals);
M_mat = X_mat(idx_j, :)' * (X_mat(idx_j, :) .* weights);

% Precision: S * Delta * M^{-1} * Delta * S
Prec = S * Delta * ((M_mat + 1e-9 * eye(d)) \ Delta) * S;

end
