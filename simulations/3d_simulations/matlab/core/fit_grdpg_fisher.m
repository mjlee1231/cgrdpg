function [X_opt, Z_opt, fval, exitflag, output, S_estimated] = fit_grdpg_fisher(A, B, d, p, tau, maxit, tol)
% FIT_GRDPG_FISHER Fit GRDPG using Fisher scoring (matching R implementation)
%
% Implements Fisher scoring with coordinate descent to match R's fit_grdpg_cov:
%   1. Initialize with augmented ASE
%   2. Iterate:
%      a. Sweep through vertices with Fisher scoring updates
%      b. Refit Z given updated X
%      c. Check convergence on max row change
%
% Inputs:
%   A     - (n x n) adjacency matrix
%   B     - (p_cov x n) covariate matrix
%   d     - embedding dimension
%   p     - number of positive signature directions
%   tau   - smoothing parameter (default: 0.001)
%   maxit - max iterations (default: 30)
%   tol   - convergence tolerance on max row change (default: 0.005)
%
% Outputs:
%   X_opt       - (n x d) optimized latent positions
%   Z_opt       - (p_cov x d) optimized covariate coefficients
%   fval        - final objective value
%   exitflag    - 1 if converged, 0 otherwise
%   output      - structure with optimization info
%   S_estimated - (d x d) estimated signature matrix

if nargin < 5 || isempty(tau), tau = 0.001; end
if nargin < 6 || isempty(maxit), maxit = 30; end
if nargin < 7 || isempty(tol), tol = 0.005; end

n = size(A, 1);
p_cov = size(B, 1);
q = d - p;

fprintf('Fitting cgrdpg with Fisher scoring...\n');
fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
fprintf('  Max iterations: %d, tolerance: %.6f\n\n', maxit, tol);

% ===== Step 1: Initialize with augmented ASE =====
fprintf('Initializing with ASE...\n');
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);  % Diagonal augmentation

% Get unsigned ASE and signature
[X_current, S_estimated] = initialize_ase(A_aug, d, p);

fprintf('  Using signature: S = diag([');
fprintf('%+d ', diag(S_estimated)');
fprintf('])\n\n');

% Initial Z
XtX = X_current' * X_current;
Z_current = (X_current \ B')';  % Solve X' * Z' = B'

% Initial objective on ORIGINAL matrix A (not augmented)
Y_current = X_current * S_estimated;
obj_current = surrogate_objective(A, X_current, Z_current, B, S_estimated, tau);

% History
history = struct();
history.max_row_change = [];
history.objective = obj_current;

fprintf('%-5s | %-15s | %-15s\n', 'Iter', 'Objective', 'Max Row Change');
fprintf('%s\n', repmat('-', 1, 50));

exitflag = 0;
converged = false;

% ===== Step 2: Fisher scoring iterations =====
for iter = 1:maxit
    % Adaptive line search parameter (match R)
    if iter <= 8
        ls_beta = 0.8;
    else
        ls_beta = 0.4;
    end

    % Fisher sweep: update all vertices with coordinate descent
    X_new = fisher_sweep_X(A, X_current, Z_current, B, S_estimated, tau, ls_beta, 1e-4, 30);

    % Refit Z given new X
    Z_new = (X_new \ B')';

    % Compute convergence metric: max row change
    row_changes = sqrt(sum((X_new - X_current).^2, 2));
    max_row_change = max(row_changes);
    history.max_row_change(iter) = max_row_change;

    % Objective
    obj_new = surrogate_objective(A, X_new, Z_new, B, S_estimated, tau);
    history.objective(iter+1) = obj_new;

    fprintf('%5d | %+15.6e | %15.6e\n', iter, obj_new, max_row_change);

    % Check convergence
    if max_row_change < tol
        converged = true;
        exitflag = 1;
        fprintf('\nConverged: max row change (%.6e) < tol (%.6e)\n', max_row_change, tol);
        break;
    end

    % Update for next iteration
    X_current = X_new;
    Z_current = Z_new;
    obj_current = obj_new;
end

if ~converged
    fprintf('\nDid not converge: max iterations (%d) reached\n', maxit);
end

% Final output
X_opt = X_current;
Z_opt = Z_current;
fval = -obj_current;  % Return negative for consistency with minimization

output = struct();
output.iterations = iter;
output.converged = converged;
output.history = history;
output.algorithm = 'Fisher scoring with coordinate descent';

fprintf('\nOptimization complete:\n');
fprintf('  Converged: %s\n', mat2str(converged));
fprintf('  Iterations: %d\n', iter);
fprintf('  Final max row change: %.6e\n', max_row_change);
fprintf('  Final objective: %.6e\n', obj_current);

end

function [X_init, S_estimated] = initialize_ase(A_aug, d, p)
    % Initialize with ASE
    n = size(A_aug, 1);
    q = d - p;

    % Eigendecomposition
    [V, D] = eig(A_aug);
    eigvals = diag(D);

    % Sort by magnitude (absolute value)
    [~, idx] = sort(abs(eigvals), 'descend');
    eigvals = eigvals(idx);
    V = V(:, idx);

    % Signature matrix (ordered: all +1's first, then -1's)
    S_estimated = diag([ones(p, 1); -ones(q, 1)]);

    % Unsigned ASE: U |Lambda|^{1/2}
    X_init = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));
end

function f = surrogate_objective(A, X, Z, B, S, tau)
    % Compute surrogate objective value
    % f = sum((A - X*Y') .* psi(X*Y') + Psi(X*Y')) - 0.5 * ||B - Z*X'||_F^2
    % where Y = X * S

    n = size(X, 1);
    Y = X * S;

    % Network part
    XY = X * Y';  % Edge probabilities (n x n)

    % Exclude diagonal (no self-loops)
    mask = ~eye(n);
    XY_vec = XY(mask);
    A_vec = A(mask);

    % Compute psi and Psi
    [psi_vals, Psi_vals, ~, ~] = psi_functions(XY_vec, tau);

    net_part = sum((A_vec - XY_vec) .* psi_vals + Psi_vals);

    % Covariate part
    resid = B - Z * X';
    cov_part = -0.5 * sum(resid(:).^2);

    f = net_part + cov_part;
end
