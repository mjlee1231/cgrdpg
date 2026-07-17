function [X_ase_unsigned, X_ase_signed, S_estimated] = fit_ase(A, d, p)
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

% Signature matrix: R's ase_grdpg constructs it as diag(c(rep(1, p), rep(-1, q)))
% This is ORDERED: all +1's first, then all -1's (not actual eigenvalue signs!)
S_estimated = diag([ones(p, 1); -ones(q, 1)]);

% Unsigned ASE: U |Λ|^{1/2}
X_ase_unsigned = V_d * diag(sqrt(abs(lambda_d)));

% Signed ASE: X_signed = X * S_estimated (matches R exactly!)
% CRITICAL: NOT using actual eigenvalue signs, but the ordered signature matrix
% R verification: all.equal(fit$X_signed, fit$X %*% fit$sign_diag) == TRUE
X_ase_signed = X_ase_unsigned * S_estimated;

end
