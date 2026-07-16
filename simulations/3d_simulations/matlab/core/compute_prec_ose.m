function Prec = compute_prec_ose(i, X_unsigned, X_signed, eps_clip)
% COMPUTE_PREC_OSE Compute precision matrix for OSE inference
%
% Matches R's compute_prec_ose function
% Edge probability: p_ji = x̃_j^T * x̂_i (signed × unsigned)
% Precision: X̃^T * diag(w) * X̃ (using signed ASE)
%
% Inputs:
%   i           - vertex index
%   X_unsigned  - (n x d) unsigned latent positions
%   X_signed    - (n x d) signed latent positions
%   eps_clip    - clipping value for probabilities
%
% Outputs:
%   Prec - (d x d) precision matrix

% All vertices except i
idx_j = setdiff(1:size(X_unsigned, 1), i);

% Edge probability: p_ji = x̃_j^T * x̂_i
p_vals = X_signed(idx_j, :) * X_unsigned(i, :)';
p_vals = max(min(p_vals, 1 - eps_clip), eps_clip);

% Weights
w = 1 ./ (p_vals .* (1 - p_vals));

% Precision using signed ASE
Prec = X_signed(idx_j, :)' * (X_signed(idx_j, :) .* w);

end
