function X_new = compute_ose_step(A, X_unsigned, X_signed, eps_clip)
% COMPUTE_OSE_STEP One-step OSE refinement for GRDPG
%
% Matches R's compute_ose_step function
% Edge probability: p_ji = x̃_j^T * x̂_i (signed × unsigned)
% Update: x_new = x_unsigned + G^{-1} * grad
%
% Inputs:
%   A          - (n x n) adjacency matrix
%   X_unsigned - (n x d) unsigned ASE
%   X_signed   - (n x d) signed ASE
%   eps_clip   - clipping value for probabilities
%
% Outputs:
%   X_new - (n x d) updated positions

n = size(A, 1);
d = size(X_unsigned, 2);
X_new = zeros(n, d);

for i = 1:n
    x_i_unsigned = X_unsigned(i, :)';  % Column vector

    % All vertices except i
    idx_j = setdiff(1:n, i);

    % Edge probability: p_ji = x̃_j^T * x̂_i
    p_i = X_signed(idx_j, :) * x_i_unsigned;
    p_i = max(min(p_i, 1 - eps_clip), eps_clip);

    % Residuals
    resid = A(i, idx_j)' - p_i;

    % Score weights
    w_score = 1 ./ (p_i .* (1 - p_i));

    % Gradient using signed ASE
    grad = X_signed(idx_j, :)' * (resid .* w_score);

    % Fisher information using signed ASE
    G = X_signed(idx_j, :)' * (X_signed(idx_j, :) .* w_score);

    % Newton-Raphson update from unsigned ASE
    X_new(i, :) = x_i_unsigned' + (G + 1e-9 * eye(d)) \ grad;
end

end
