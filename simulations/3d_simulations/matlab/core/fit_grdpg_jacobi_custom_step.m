function [X_opt, Z_opt, fval, exitflag, output, S_estimated] = fit_grdpg_jacobi_custom_step(A, B, d, p, tau, maxit, tol, step_size_init)
% FIT_GRDPG_JACOBI_CUSTOM_STEP Jacobi coordinate ascent with custom step size
%
% Same as fit_grdpg_jacobi but allows testing different step sizes
%
% Inputs:
%   step_size_init - initial step size to use (will reduce adaptively)

if nargin < 5 || isempty(tau), tau = 0.001; end
if nargin < 6 || isempty(maxit), maxit = 30; end
if nargin < 7 || isempty(tol), tol = 0.005; end
if nargin < 8 || isempty(step_size_init), step_size_init = 0.1; end

n = size(A, 1);
p_cov = size(B, 1);
q = d - p;

% Suppress output for speed during testing
verbose = false;

if verbose
    fprintf('Fitting cgrdpg with Jacobi (step_init=%.3f)...\n', step_size_init);
    fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
end

% Initialize with ASE
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

[X_current, S_estimated] = initialize_ase(A_aug, d, p);

% History
history = struct();
history.max_row_change = [];
history.objective = [];
history.step_size = [];

exitflag = 0;
converged = false;

% Jacobi coordinate ascent iterations
for iter = 1:maxit
    % Current Y and Z
    Y_current = X_current * S_estimated;
    Z_current = (X_current \ B')';

    % Adaptive step size schedule based on initial value
    % Reduce step size over iterations
    if iter <= 5
        step_size = step_size_init;
    elseif iter <= 15
        step_size = step_size_init * 0.5;  % Half after iter 5
    else
        step_size = step_size_init * 0.1;  % 1/10 after iter 15
    end

    % Vectorized update
    X_new = jacobi_sweep_vectorized(A, X_current, Y_current, Z_current, B, S_estimated, tau, step_size);

    % Refit Z
    Z_new = (X_new \ B')';

    % Compute objective
    fval = compute_objective(A, X_new, Y_current, Z_current, B, tau);
    history.objective(iter) = fval;

    % Convergence metric
    row_changes = sqrt(sum((X_new - X_current).^2, 2));
    max_row_change = max(row_changes);
    history.max_row_change(iter) = max_row_change;
    history.step_size(iter) = step_size;

    % Check convergence
    if max_row_change < tol
        converged = true;
        exitflag = 1;
        break;
    end

    % Check for numerical explosion
    if ~isfinite(fval) || max_row_change > 1e6
        converged = false;
        exitflag = -1;
        break;
    end

    % Update
    X_current = X_new;
    Z_current = Z_new;
end

X_opt = X_current;
Z_opt = Z_current;

output = struct();
output.iterations = iter;
output.converged = converged;
output.history = history;
output.algorithm = sprintf('Jacobi (step_init=%.3f)', step_size_init);

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
    n = size(X, 1);
    d = size(X, 2);

    Y_t = Y';
    ZtZ = Z' * Z;
    X_new = zeros(n, d);

    for i = 1:n
        x_i = X(i, :)';

        % Edge probabilities
        s = Y * x_i;
        w = dpsi(s, tau);
        w(i) = 0;

        % Network score
        r = (A(i, :)' - s) .* w;
        s_net = Y_t * r;

        % Covariate score
        b_i = B(:, i);
        resid = b_i - Z * x_i;
        s_cov = Z' * resid;

        S_score = s_net + s_cov;

        % Fisher information
        s_clipped = max(min(s, 1 - tau), tau);
        w_fisher = dpsi(s_clipped, tau);
        w_fisher(i) = 0;

        G_net = Y_t * (Y .* w_fisher);
        G = G_net + ZtZ;

        % Newton direction
        [R, flag] = chol(G);
        if flag == 0
            p = R \ (R' \ S_score);
        else
            p = (G + 1e-9 * eye(d)) \ S_score;
        end

        % Check if ascent
        sp = S_score' * p;
        if ~isfinite(sp) || sp <= 0
            p = S_score;
            sp = S_score' * p;
            if ~isfinite(sp) || sp <= 0
                X_new(i, :) = x_i';
                continue;
            end
        end

        % Update with custom step size
        X_new(i, :) = (x_i + step_size * p)';
    end
end

function obj = compute_objective(A, X, Y, Z, B, tau)
    n = size(X, 1);
    XY = X * Y';

    mask = ~eye(n);
    XY_vec = XY(mask);
    A_vec = A(mask);

    [psi_vals, Psi_vals, ~, ~] = psi_functions(XY_vec, tau);
    net_part = sum((A_vec - XY_vec) .* psi_vals + Psi_vals);

    resid = B - Z * X';
    cov_part = -0.5 * sum(resid(:).^2);

    obj = net_part + cov_part;
end
