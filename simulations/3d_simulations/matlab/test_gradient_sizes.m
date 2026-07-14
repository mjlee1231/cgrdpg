% Test gradient correctness at different problem sizes
clear; clc;

% Add core folder to path
addpath('core');

fprintf('Testing Gradient at Multiple Problem Sizes\n');
fprintf('========================================\n\n');

% Test different sizes
test_sizes = [5, 10, 20, 30, 50];
tau = 0.001;

results = struct();

for size_idx = 1:length(test_sizes)
    n = test_sizes(size_idx);
    p_cov = max(5, round(n/2));
    d = 2;

    fprintf('========================================\n');
    fprintf('Testing n=%d, p_cov=%d, d=%d\n', n, p_cov, d);
    fprintf('========================================\n\n');

    rng(42);  % Same seed for consistency

    % Generate simple data (2D semicircle)
    theta = pi * (1:n)' / (n - 1);
    X_true = [0.28 * sin(theta) + 0.42, ...
              0.28 * cos(theta) + 0.42];
    Z_true = randn(p_cov, d) * 0.1;
    S = eye(d);  % Identity (all positive)
    Y_true = X_true * S;
    P_net = X_true * Y_true';

    % Print probability range
    fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_net(:)), max(P_net(:)));

    A = double(rand(n) < P_net);
    A = triu(A, 1);
    A = A + A';
    A(1:n+1:end) = 0;  % Ensure diagonal is 0
    B = Z_true * X_true' + randn(p_cov, n) * 0.1;

    % Initialize with ASE
    A_aug = A;
    A_aug(1:n+1:end) = sum(A, 2) / (n - 1);
    [V, D] = eig(A_aug);
    eigvals = diag(D);
    [~, idx] = sort(abs(eigvals), 'descend');
    eigvals_sorted = eigvals(idx);
    X_init = V(:, idx(1:d)) * diag(sqrt(abs(eigvals_sorted(1:d))));

    % Estimate signature from ASE eigenvalues
    S_estimated = diag(sign(eigvals_sorted(1:d)));

    % Fix Y_hat and Z_hat
    Y_hat = X_init * S_estimated;
    Z_hat = (X_init \ B')';

    % Pack X only
    x0 = X_init(:);

    % Define objective
    objective = @(x) surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau);

    % Compute analytical gradient
    [f, g_analytical] = objective(x0);

    % Compute finite difference gradient (sample 10 components randomly)
    eps_fd = 1e-7;
    num_samples = min(10, length(x0));
    sample_idx = randperm(length(x0), num_samples);

    ratios = zeros(num_samples, 1);
    abs_diffs = zeros(num_samples, 1);

    for i = 1:num_samples
        idx = sample_idx(i);
        x_plus = x0;
        x_plus(idx) = x_plus(idx) + eps_fd;
        f_plus = objective(x_plus);

        x_minus = x0;
        x_minus(idx) = x_minus(idx) - eps_fd;
        f_minus = objective(x_minus);

        g_fd = (f_plus - f_minus) / (2 * eps_fd);

        ratios(i) = g_analytical(idx) / g_fd;
        abs_diffs(i) = abs(g_analytical(idx) - g_fd);
    end

    % Store results
    results(size_idx).n = n;
    results(size_idx).mean_ratio = mean(ratios);
    results(size_idx).std_ratio = std(ratios);
    results(size_idx).min_ratio = min(ratios);
    results(size_idx).max_ratio = max(ratios);
    results(size_idx).max_abs_diff = max(abs_diffs);
    results(size_idx).mean_abs_diff = mean(abs_diffs);

    fprintf('Sample gradient check (%d random components):\n', num_samples);
    fprintf('  Mean ratio (analytical/fd): %.6f\n', mean(ratios));
    fprintf('  Std ratio: %.6f\n', std(ratios));
    fprintf('  Min ratio: %.6f\n', min(ratios));
    fprintf('  Max ratio: %.6f\n', max(ratios));
    fprintf('  Max abs diff: %.6e\n', max(abs_diffs));
    fprintf('  Mean abs diff: %.6e\n\n', mean(abs_diffs));

    if abs(mean(ratios) - 1.0) < 0.01
        fprintf('✓ Gradient appears CORRECT at n=%d\n\n', n);
    else
        fprintf('✗ Gradient has errors at n=%d (ratio = %.3f)\n\n', n, mean(ratios));
    end
end

%% Summary
fprintf('========================================\n');
fprintf('Summary Across All Sizes\n');
fprintf('========================================\n\n');
fprintf('%-6s | %-12s | %-12s | %-12s\n', 'n', 'Mean Ratio', 'Std Ratio', 'Max Abs Diff');
fprintf('%s\n', repmat('-', 1, 50));
for i = 1:length(results)
    fprintf('%-6d | %12.6f | %12.6f | %12.6e\n', ...
        results(i).n, results(i).mean_ratio, results(i).std_ratio, results(i).max_abs_diff);
end
fprintf('\n');

%% Nested function
function [f, g] = surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau)
    X = reshape(x, n, d);
    S_mat = X * Y_hat';
    is_off_diag = ~eye(n);

    [psi_val, Psi_val, dpsi_val] = psi_functions(S_mat, tau);

    net_obj = sum((A(is_off_diag) - S_mat(is_off_diag)) .* psi_val(is_off_diag) + Psi_val(is_off_diag));

    B_pred = Z_hat * X';
    cov_obj = -0.5 * sum((B - B_pred).^2, 'all');

    f = -(net_obj + cov_obj);

    if nargout > 1
        W_core = (A - S_mat) .* dpsi_val;
        W_core(~is_off_diag) = 0;
        W_net = W_core;
        grad_X_net = -W_net * Y_hat;
        grad_X_cov = -(B - B_pred)' * Z_hat;
        grad_X = grad_X_net + grad_X_cov;
        g = grad_X(:);
    end
end
