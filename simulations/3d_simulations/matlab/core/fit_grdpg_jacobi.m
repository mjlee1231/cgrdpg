function [X_opt, Z_opt, fval, exitflag, output, S_estimated] = fit_grdpg_jacobi(A, B, d, p, tau, maxit, tol)
% FIT_GRDPG_JACOBI Vectorized Jacobi block coordinate ascent for GRDPG
%
% Implements vectorized Jacobi-style coordinate descent:
% - Computes updates for ALL vertices based on CURRENT state (vectorized)
% - Applies all updates simultaneously (Jacobi instead of Gauss-Seidel)
% - Middle ground between R's sequential updates and MATLAB's batch surrogate
%
% Inputs:
%   A     - (n x n) adjacency matrix
%   B     - (p_cov x n) covariate matrix
%   d     - embedding dimension
%   p     - number of positive eigenvalues
%   tau   - smoothing parameter (default: 0.001)
%   maxit - max iterations (default: 30)
%   tol   - convergence tolerance on max row change (default: 0.005)
%
% Outputs:
%   X_opt       - (n x d) optimized latent positions
%   Z_opt       - (p_cov x d) covariate coefficients
%   fval        - final objective value
%   exitflag    - 1 if converged, 0 otherwise
%   output      - optimization info
%   S_estimated - (d x d) signature matrix

if nargin < 5 || isempty(tau), tau = 0.001; end
if nargin < 6 || isempty(maxit), maxit = 30; end
if nargin < 7 || isempty(tol), tol = 0.005; end

n = size(A, 1);
p_cov = size(B, 1);
q = d - p;

fprintf('Fitting cgrdpg with vectorized Jacobi coordinate ascent...\n');
fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
fprintf('  Max iterations: %d, tolerance: %.6f\n\n', maxit, tol);

% Initialize with ASE
fprintf('Initializing with ASE...\n');
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

[X_current, S_estimated] = initialize_ase(A_aug, d, p);

fprintf('  Using signature: S = diag([');
fprintf('%+d ', diag(S_estimated)');
fprintf('])\n\n');

% History
history = struct();
history.max_row_change = [];
history.objective = [];
history.step_size = [];

fprintf('%-5s | %-15s | %-15s | %-10s\n', 'Iter', 'Objective', 'Max Row Change', 'Step Size');
fprintf('%s\n', repmat('-', 1, 65));

exitflag = 0;
converged = false;

% Jacobi coordinate ascent iterations
for iter = 1:maxit
    % Current Y and Z
    Y_current = X_current * S_estimated;
    Z_current = (X_current \ B')';

    % Adaptive step size: very conservative to handle extreme edge probs
    % When edge probs are near 1 (like 0.9971), Fisher info explodes
    % Use very small steps to prevent numerical catastrophe
    if iter <= 5
        step_size = 0.01;  % Initial phase: very small step
    elseif iter <= 15
        step_size = 0.005;  % Middle phase: even smaller
    else
        step_size = 0.001;  % Final phase: tiny refinement steps
    end

    % Vectorized update: compute Newton direction for ALL vertices
    X_new = jacobi_sweep_vectorized(A, X_current, Y_current, Z_current, B, S_estimated, tau, step_size);

    % Refit Z
    Z_new = (X_new \ B')';

    % Compute objective (for monitoring)
    fval = compute_objective(A, X_new, Y_current, Z_current, B, tau);
    history.objective(iter) = fval;

    % Convergence metric
    row_changes = sqrt(sum((X_new - X_current).^2, 2));
    max_row_change = max(row_changes);
    history.max_row_change(iter) = max_row_change;
    history.step_size(iter) = step_size;

    fprintf('%5d | %+15.6e | %15.6e | %10.3f\n', iter, fval, max_row_change, step_size);

    % Check convergence
    if max_row_change < tol
        converged = true;
        exitflag = 1;
        fprintf('\nConverged: max row change (%.6e) < tol (%.6e)\n', max_row_change, tol);
        break;
    end

    % Update
    X_current = X_new;
    Z_current = Z_new;
end

if ~converged
    fprintf('\nDid not converge: max iterations (%d) reached\n', maxit);
end

X_opt = X_current;
Z_opt = Z_current;

output = struct();
output.iterations = iter;
output.converged = converged;
output.history = history;
output.algorithm = 'Vectorized Jacobi coordinate ascent';

fprintf('\nOptimization complete:\n');
fprintf('  Converged: %s\n', mat2str(converged));
fprintf('  Iterations: %d\n', iter);
fprintf('  Final max row change: %.6e\n', max_row_change);
fprintf('  Final objective: %.6e\n', fval);

end

function [X_init, S_estimated] = initialize_ase(A_aug, d, p)
    n = size(A_aug, 1);
    q = d - p;

    [V, D] = eig(A_aug);
    eigvals = diag(D);

    [~, idx] = sort(abs(eigvals), 'descend');
    eigvals = eigvals(idx);
    V = V(:, idx);

    S_estimated = diag([ones(p, 1); -ones(q, 1)]);
    X_init = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));
end

function X_new = jacobi_sweep_vectorized(A, X, Y, Z, B, S, tau, step_size)
% Vectorized Jacobi sweep: compute updates for ALL vertices simultaneously
% Based on current state (X, Y, Z), compute Newton directions for all vertices
% Then apply all updates at once (Jacobi-style, not Gauss-Seidel)
%
% Inputs:
%   step_size - damping factor for Newton step (0 < step_size <= 1)

n = size(X, 1);
d = size(X, 2);

% Pre-compute shared quantities
Y_t = Y';  % (d x n)
ZtZ = Z' * Z;  % (d x d)

% Allocate output
X_new = zeros(n, d);

% Vectorized computation for all vertices
for i = 1:n
    x_i = X(i, :)';  % (d x 1)

    % Edge probabilities: s_j = x_i^T * y_j for all j != i
    s = Y * x_i;  % (n x 1)
    w = dpsi(s, tau);
    w(i) = 0;  % Exclude self

    % Network score
    r = (A(i, :)' - s) .* w;
    s_net = Y_t * r;  % (d x 1)

    % Covariate score
    b_i = B(:, i);
    resid = b_i - Z * x_i;
    s_cov = Z' * resid;  % (d x 1)

    S_score = s_net + s_cov;  % Total score (d x 1)

    % Fisher information (clip for stability)
    s_clipped = max(min(s, 1 - tau), tau);
    w_fisher = dpsi(s_clipped, tau);
    w_fisher(i) = 0;

    G_net = Y_t * (Y .* w_fisher);  % (d x d)
    G = G_net + ZtZ;

    % Newton direction
    [R, flag] = chol(G);
    if flag == 0
        p = R \ (R' \ S_score);
    else
        p = (G + 1e-9 * eye(d)) \ S_score;
    end

    % Check if ascent direction
    sp = S_score' * p;
    if ~isfinite(sp) || sp <= 0
        p = S_score;
        sp = S_score' * p;
        if ~isfinite(sp) || sp <= 0
            % Skip update for this vertex
            X_new(i, :) = x_i';
            continue;
        end
    end

    % Jacobi update with adaptive damping
    % step_size is passed from outer loop and decreases over iterations
    X_new(i, :) = (x_i + step_size * p)';
end

end

function obj = compute_objective(A, X, Y, Z, B, tau)
% Compute surrogate objective value
n = size(X, 1);
XY = X * Y';  % Edge probabilities

% Network part (exclude diagonal)
mask = ~eye(n);
XY_vec = XY(mask);
A_vec = A(mask);

[psi_vals, Psi_vals, ~, ~] = psi_functions(XY_vec, tau);
net_part = sum((A_vec - XY_vec) .* psi_vals + Psi_vals);

% Covariate part
resid = B - Z * X';
cov_part = -0.5 * sum(resid(:).^2);

obj = net_part + cov_part;
end
