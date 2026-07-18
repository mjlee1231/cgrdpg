% Fine-tune step size in narrow range 0.015 - 0.025
% Previous results: 0.020 worked (SSE=59.39), 0.030 exploded

clear; clc;
addpath('core');

fprintf('Fine-Tuning Jacobi Step Size (0.015 - 0.025)\n');
fprintf('=============================================\n\n');

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

% Fine-grained step sizes around 0.020
step_sizes_test = [0.012, 0.015, 0.017, 0.019, 0.020, 0.021, 0.022, 0.023, 0.024, 0.025, 0.027];

results = struct();
results.step_size = [];
results.sse = [];
results.time = [];
results.converged = [];
results.iterations = [];
results.status = {};

fprintf('Testing %d step sizes in range [0.012, 0.027]...\n\n', length(step_sizes_test));
fprintf('%-10s | %-12s | %-10s | %-8s | %-15s\n', 'Step Size', 'SSE', 'Time (s)', 'Iters', 'Status');
fprintf('%s\n', repmat('-', 1, 75));

for idx = 1:length(step_sizes_test)
    step_init = step_sizes_test(idx);

    try
        t0 = tic;
        [X_jac, ~, ~, exitflag, output, ~] = ...
            fit_grdpg_jacobi_custom_step(A, B, d, p, tau, maxit, tol, step_init);
        time_elapsed = toc(t0);

        [X_aligned, ~] = procrustes_align(X_jac, X0);
        sse = sum((X_aligned - X0).^2, 'all');

        % Classify result
        if ~isfinite(sse) || sse > 1e6
            status = 'EXPLODED';
        elseif sse < 51.15
            status = '✓ BETTER';
        elseif sse < 60
            status = 'CLOSE';
        else
            status = 'WORSE';
        end

        results.step_size(idx) = step_init;
        results.sse(idx) = sse;
        results.time(idx) = time_elapsed;
        results.converged(idx) = exitflag;
        results.iterations(idx) = output.iterations;
        results.status{idx} = status;

        fprintf('%-10.3f | %-12.2f | %-10.1f | %-8d | %-15s\n', ...
            step_init, sse, time_elapsed, output.iterations, status);

    catch ME
        fprintf('%-10.3f | %-12s | %-10s | %-8s | %-15s\n', ...
            step_init, 'ERROR', '-', '-', 'ERROR');

        results.step_size(idx) = step_init;
        results.sse(idx) = NaN;
        results.time(idx) = NaN;
        results.converged(idx) = false;
        results.iterations(idx) = 0;
        results.status{idx} = 'ERROR';
    end
end

fprintf('%s\n', repmat('-', 1, 75));

%% Analysis
valid_idx = isfinite(results.sse) & (results.sse < 1e6);
if any(valid_idx)
    valid_sse = results.sse(valid_idx);
    valid_steps = results.step_size(valid_idx);

    [best_sse, best_local_idx] = min(valid_sse);
    best_step = valid_steps(best_local_idx);

    fprintf('\n========================================\n');
    fprintf('Results Summary\n');
    fprintf('========================================\n');
    fprintf('Best step size: %.3f → SSE = %.2f\n', best_step, best_sse);

    fprintf('\nAll stable results:\n');
    [sorted_sse, sort_idx] = sort(valid_sse);
    for i = 1:min(5, length(sorted_sse))
        fprintf('  %.3f → SSE = %.2f\n', valid_steps(sort_idx(i)), sorted_sse(i));
    end

    fprintf('\n========================================\n');
    fprintf('Comparison with Baselines\n');
    fprintf('========================================\n');
    fprintf('R Fisher scoring:            SSE = 8.23\n');
    fprintf('Batch Trust-Region (MATLAB): SSE = 51.15\n');
    fprintf('Best Jacobi (step=%.3f):    SSE = %.2f\n', best_step, best_sse);

    if best_sse < 51.15
        improvement = 100 * (51.15 - best_sse) / 51.15;
        fprintf('\n✓✓✓ SUCCESS! %.1f%% improvement over batch!\n', improvement);
        fprintf('Recommended: Use Jacobi with step_size = %.3f\n', best_step);
    elseif best_sse < 55
        degradation = 100 * (best_sse - 51.15) / 51.15;
        fprintf('\n~ Close but not better (%.1f%% worse than batch)\n', degradation);
        fprintf('Recommendation: Stick with batch surrogate method\n');
    else
        degradation = 100 * (best_sse - 51.15) / 51.15;
        fprintf('\n✗ Still significantly worse (%.1f%% worse than batch)\n', degradation);
        fprintf('Recommendation: Abandon Jacobi, use batch surrogate\n');
    end

    % Plot if possible
    if length(valid_sse) > 3
        fprintf('\nSSE vs Step Size trend:\n');
        fprintf('  Smallest step tested: %.3f → SSE = %.2f\n', ...
            min(valid_steps), valid_sse(valid_steps == min(valid_steps)));
        fprintf('  Largest step tested:  %.3f → SSE = %.2f\n', ...
            max(valid_steps), valid_sse(valid_steps == max(valid_steps)));
    end
else
    fprintf('\n✗✗✗ All step sizes failed!\n');
    fprintf('Recommendation: Abandon Jacobi method entirely\n');
end
