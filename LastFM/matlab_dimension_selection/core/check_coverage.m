function covered = check_coverage(err, Prec, chi2_crit, scale)
% CHECK_COVERAGE Check if true position is within confidence region
%
% Tests whether the true latent position falls within the (1-alpha)
% confidence ellipsoid around the estimated position using the
% chi-squared statistic.
%
% Inputs:
%   err       - (1 x d) error vector: X_true[i,:] - X_est[i,:]
%   Prec      - (d x d) precision matrix (Fisher information)
%   chi2_crit - chi-squared critical value for desired confidence level
%   scale     - optional scaling factor (default: 1.0)
%
% Outputs:
%   covered - logical, true if within confidence region, false otherwise
%             returns NaN if precision matrix is not positive definite
%
% Test statistic:
%   Q = scale * err * Prec * err'
%   covered = (Q <= chi2_crit)

if nargin < 4
    scale = 1.0;
end

% Check if precision matrix is positive definite
eigvals = eig(Prec);
if min(eigvals) < 1e-10
    covered = NaN;
    return;
end

% Compute chi-squared test statistic
chi2_stat = scale * (err * Prec * err');

% Check coverage
covered = (chi2_stat <= chi2_crit);

end
