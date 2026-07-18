% Test different step sizes for Jacobi method to find optimal
% Goal: Find step size that balances numerical stability and optimization quality

clear; clc;
addpath('core');

fprintf('Testing Jacobi Method with Different Step Sizes\n');
fprintf('===============================================\n\n');

% Parameters (rep 1)
n = 1000;
p_cov = 500;
d = 3;
p = 2;
tau = 0.001;
maxit = 30;
tol = 0.005;

% Set seed
rng(598 + 1);

% Generate data
t = (1:n)' / n;
X0 = [0.15 * sin(2*pi*t) + 0.6, ...
      0.15 * cos(2*pi*t) + 0.6, ...
      0.15 * cos(4*pi*t)];

S = diag([1, 1, -1]);
Y0 = X0 * S;
Z0 = randn(p_cov, d);
P = X0 * Y0';

A = double(rand(n) < P);
A = triu(A, 1);
A = A + A';
A(1:n+1:end) = 0;
B = Z0 * X0' + randn(p_cov, n);

fprintf('Data: n=%d, p_cov=%d, d=%d\n', n, p_cov, d);
fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P(:)), max(P(:)));

% Test different initial step sizes
% Strategy: Start with different values, keep adaptive reduction
step_sizes_initial = [0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25];

results = struct();
results.step_size = [];
results.sse = [];
results.time = [];
results.converged = [];
results.iterations = [];
results.status = {};

fprintf('Testing %d different initial step sizes...\n\n', length(step_sizes_initial));
fprintf('%-10s | %-10s | %-10s | %-10s | %-15s\n', 'Step Size', 'SSE', 'Time (s)', 'Iters', 'Status');
fprintf('%s\n', repmat('-', 1, 70));

for idx = 1:length(step_sizes_initial)
    step_init = step_sizes_initial(idx);

    try
        % Modify fit_grdpg_jacobi to accept custom step size schedule
        % For now, we'll test with constant step size
        t0 = tic;
        [X_jac, ~, ~, exitflag, output, ~] = ...
            fit_grdpg_jacobi_custom_step(A, B, d, p, tau, maxit, tol, step_init);
        time_elapsed = toc(t0);

        [X_aligned, ~] = procrustes_align(X_jac, X0);
        sse = sum((X_aligned - X0).^2, 'all');

        % Check for numerical issues
        if ~isfinite(sse) || sse > 1e6
            status = 'EXPLODED';
        elseif sse < 51.15
            status = 'BETTER';
        else
            status = 'WORSE';
        end

        results.step_size(idx) = step_init;
        results.sse(idx) = sse;
        results.time(idx) = time_elapsed;
        results.converged(idx) = exitflag;
        results.iterations(idx) = output.iterations;
        results.status{idx} = status;

        fprintf('%-10.3f | %-10.2f | %-10.1f | %-10d | %-15s\n', ...
            step_init, sse, time_elapsed, output.iterations, status);

    catch ME
        fprintf('%-10.3f | %-10s | %-10s | %-10s | %-15s\n', ...
            step_init, 'ERROR', '-', '-', ME.message(1:min(15,end)));

        results.step_size(idx) = step_init;
        results.sse(idx) = NaN;
        results.time(idx) = NaN;
        results.converged(idx) = false;
        results.iterations(idx) = 0;
        results.status{idx} = 'ERROR';
    end
end

fprintf('%s\n', repmat('-', 1, 70));

%% Find best result
valid_idx = isfinite(results.sse) & (results.sse < 1e6);
if any(valid_idx)
    valid_sse = results.sse(valid_idx);
    valid_steps = results.step_size(valid_idx);

    [best_sse, best_local_idx] = min(valid_sse);
    best_step = valid_steps(best_local_idx);

    fprintf('\n========================================\n');
    fprintf('Best Result\n');
    fprintf('========================================\n');
    fprintf('Step size: %.3f\n', best_step);
    fprintf('SSE:       %.2f\n', best_sse);

    fprintf('\nComparison:\n');
    fprintf('  R Fisher scoring:         SSE = 8.23\n');
    fprintf('  Batch Trust-Region:       SSE = 51.15\n');
    fprintf('  Best Jacobi (step=%.3f): SSE = %.2f\n', best_step, best_sse);

    if best_sse < 51.15
        improvement = 100 * (51.15 - best_sse) / 51.15;
        fprintf('\n✓ SUCCESS: %.1f%% improvement over batch!\n', improvement);
    else
        degradation = 100 * (best_sse - 51.15) / 51.15;
        fprintf('\n✗ Still worse than batch by %.1f%%\n', degradation);
    end
else
    fprintf('\n✗ All step sizes failed!\n');
end
