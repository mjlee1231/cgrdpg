function [X_ase_unsigned, X_ase_signed, S_estimated, lambda_d, V_d] = fit_ase(A, d, p)
% FIT_ASE Adjacency Spectral Embedding for GRDPG
%
% Matches R's ase_grdpg function
%
% Inputs:
%   A - (n x n) adjacency matrix
%   d - embedding dimension
%   p - number of positive eigenvalues for signature construction
%
% Outputs:
%   X_ase_unsigned - (n x d) unsigned ASE: U |Λ|^{1/2}
%   X_ase_signed   - (n x d) signed ASE: U |Λ|^{1/2} sign(Λ)
%   S_estimated    - (d x d) signature matrix from actual eigenvalue signs
%   lambda_d       - (d x 1) top d eigenvalues (sorted by magnitude)
%   V_d            - (n x d) top d eigenvectors

n = size(A, 1);
q = d - p;

% Augmented adjacency matrix
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

% Eigendecomposition
[V, D] = eig(A_aug);
eigvals = diag(D);

% Sort by magnitude (absolute value) - matches R's eigs_sym(..., which="LM")
[~, idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(idx);
V = V(:, idx);

% Take top d eigenvectors
V_d = V(:, 1:d);
lambda_d = eigvals(1:d);

% Use ACTUAL eigenvalue signs (matches R's implementation)
% R: X_signed <- U %*% diag(sqrt(abs(vals)) * sign(vals))
S_estimated = diag(sign(lambda_d));

% Unsigned ASE: U |Λ|^{1/2}
X_ase_unsigned = V_d * diag(sqrt(abs(lambda_d)));

% Signed ASE: X_signed = U |Λ|^{1/2} sign(Λ)
% Uses actual observed eigenvalue signs, not ordered signature
X_ase_signed = X_ase_unsigned * S_estimated;

end
