function [X_opt, Z_opt, fval, exitflag, output, S_estimated] = fit_grdpg_fisher_fast(A, B, d, p, tau, maxit, tol)
% FIT_GRDPG_FISHER_FAST Fast Fisher scoring without expensive line search
%
% Uses fixed step size instead of backtracking line search for speed
% Should run in ~2-3 minutes instead of hours

if nargin < 5 || isempty(tau), tau = 0.001; end
if nargin < 6 || isempty(maxit), maxit = 30; end
if nargin < 7 || isempty(tol), tol = 0.005; end

n = size(A, 1);
p_cov = size(B, 1);
q = d - p;

fprintf('Fitting cgrdpg with fast Fisher scoring (no line search)...\n');
fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
fprintf('  Max iterations: %d, tolerance: %.6f\n\n', maxit, tol);

% Initialize with augmented ASE
fprintf('Initializing with ASE...\n');
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

[X_current, S_estimated] = initialize_ase(A_aug, d, p);

fprintf('  Using signature: S = diag([');
fprintf('%+d ', diag(S_estimated)');
fprintf('])\n\n');

% Initial Z
Z_current = (X_current \ B')';

% History
history = struct();
history.max_row_change = [];

fprintf('%-5s | %-15s\n', 'Iter', 'Max Row Change');
fprintf('%s\n', repmat('-', 1, 30));

exitflag = 0;
converged = false;

% Fisher scoring iterations
for iter = 1:maxit
    % Adaptive step size (match R's adaptive ls_beta)
    if iter <= 8
        step_size = 0.8;
    else
        step_size = 0.4;
    end

    % Fisher sweep with fixed step size (no line search)
    X_new = fisher_sweep_X_fast(A, X_current, Z_current, B, S_estimated, tau, step_size);

    % Refit Z
    Z_new = (X_new \ B')';

    % Convergence metric
    row_changes = sqrt(sum((X_new - X_current).^2, 2));
    max_row_change = max(row_changes);
    history.max_row_change(iter) = max_row_change;

    fprintf('%5d | %15.6e\n', iter, max_row_change);

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
fval = 0;  % Don't compute expensive objective

output = struct();
output.iterations = iter;
output.converged = converged;
output.history = history;
output.algorithm = 'Fast Fisher scoring (fixed step size)';

fprintf('\nOptimization complete:\n');
fprintf('  Converged: %s\n', mat2str(converged));
fprintf('  Iterations: %d\n', iter);
fprintf('  Final max row change: %.6e\n', max_row_change);

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

function X_new = fisher_sweep_X_fast(A, X, Z, B, S, tau, step_size)
% Fast Fisher sweep with FIXED step size (no line search)

n = size(X, 1);
d = size(X, 2);

X_cur = X;
Y_cur = X_cur * S;
Y_cur_t = Y_cur';
ZtZ = Z' * Z;

for i = 1:n
    x_i = X_cur(i, :)';

    % Gradient computation
    s = Y_cur * x_i;
    w = dpsi(s, tau);
    w(i) = 0;

    r = (A(i, :)' - s) .* w;
    s_net = Y_cur_t * r;

    b_i = B(:, i);
    resid = b_i - Z * x_i;
    s_cov = Z' * resid;

    S_score = s_net + s_cov;

    % Fisher information
    s_clipped = max(min(s, 1 - tau), tau);
    w_fisher = dpsi(s_clipped, tau);
    w_fisher(i) = 0;

    G_net = Y_cur_t * (Y_cur .* w_fisher);
    G = G_net + ZtZ;

    % Solve for Newton direction
    [R, flag] = chol(G);
    if flag == 0
        p = R \ (R' \ S_score);
    else
        p = G \ S_score;
    end

    % Check if ascent direction
    sp = S_score' * p;
    if ~isfinite(sp) || sp <= 0
        p = S_score;
        sp = S_score' * p;
        if ~isfinite(sp) || sp <= 0
            continue;
        end
    end

    % FIXED step size (no expensive line search!)
    X_cur(i, :) = (x_i + step_size * p)';
    Y_cur(i, :) = X_cur(i, :) * S;
    Y_cur_t = Y_cur';
end

X_new = X_cur;
end
