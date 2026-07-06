% Test gradient correctness using finite differences
clear; clc;

fprintf('Testing Gradient Correctness\n');
fprintf('========================================\n\n');

%% Small test problem
n = 50;
p_cov = 25;
d = 3;
p = 2;
tau = 0.001;
rng(42);

% Generate simple data
t = (1:n)' / n;
X_true = [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)];
Z_true = randn(p_cov, d) * 0.1;
S = diag([1, 1, -1]);
Y_true = X_true * S;
P_net = X_true * Y_true';
A = double(rand(n) < P_net);
A = triu(A, 1);
A = A + A';
B = Z_true * X_true' + randn(p_cov, n) * 0.1;

% Initialize with ASE
A_aug = A;
A_aug(1:n+1:end) = sum(A, 2) / (n - 1);
[V, D] = eig(A_aug);
eigvals = diag(D);
[~, idx] = sort(abs(eigvals), 'descend');
X_init = V(:, idx(1:d)) * diag(sqrt(abs(eigvals(idx(1:d)))));
Z_init = B * X_init / (X_init' * X_init);

% Pack parameters
x0 = [X_init(:); Z_init(:)];

% Define objective
objective = @(x) surrogate_objective_gradient(x, A, B, S, n, d, p_cov, tau);

%% Test gradient with MATLAB's built-in checker
fprintf('Testing gradient with finite differences...\n\n');

options = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'CheckGradients', true, ...
    'FiniteDifferenceType', 'central', ...
    'Display', 'off');

% Run one iteration to trigger gradient check
fprintf('Running fminunc with gradient checking enabled...\n');
[~, ~, ~, output] = fminunc(objective, x0, options);

fprintf('\nGradient check complete.\n');
fprintf('If gradients are correct, you should see small differences above.\n');
fprintf('If differences are large (>1e-6), there is a bug in the gradient!\n');

%% Also compute gradient manually
fprintf('\n========================================\n');
fprintf('Manual gradient check at initial point\n');
fprintf('========================================\n\n');

[f, g_analytical] = objective(x0);

% Finite difference gradient
eps_fd = 1e-7;
g_fd = zeros(size(g_analytical));
for i = 1:length(x0)
    x_plus = x0;
    x_plus(i) = x_plus(i) + eps_fd;
    f_plus = objective(x_plus);

    x_minus = x0;
    x_minus(i) = x_minus(i) - eps_fd;
    f_minus = objective(x_minus);

    g_fd(i) = (f_plus - f_minus) / (2 * eps_fd);
end

% Compare
diff = abs(g_analytical - g_fd);
rel_diff = diff ./ (abs(g_fd) + 1e-10);

fprintf('Objective value: %.6e\n', f);
fprintf('Gradient norm (analytical): %.6e\n', norm(g_analytical));
fprintf('Gradient norm (finite diff): %.6e\n', norm(g_fd));
fprintf('\n');
fprintf('Max absolute difference: %.6e\n', max(diff));
fprintf('Max relative difference: %.6e\n', max(rel_diff));
fprintf('Mean relative difference: %.6e\n', mean(rel_diff));
fprintf('\n');

if max(rel_diff) < 1e-5
    fprintf('✓ Gradient appears CORRECT (max rel diff < 1e-5)\n');
elseif max(rel_diff) < 1e-3
    fprintf('⚠ Gradient may have minor errors (max rel diff < 1e-3)\n');
else
    fprintf('✗ Gradient has MAJOR ERRORS (max rel diff > 1e-3)!\n');
    fprintf('\nShowing worst 5 components:\n');
    [~, worst_idx] = sort(rel_diff, 'descend');
    for k = 1:min(5, length(worst_idx))
        idx = worst_idx(k);
        fprintf('  Component %d: analytical=%.6e, fd=%.6e, rel_diff=%.6e\n', ...
            idx, g_analytical(idx), g_fd(idx), rel_diff(idx));
    end
end

%% Nested function matching fit_grdpg_fminunc_surrogate.m
function [f, g] = surrogate_objective_gradient(x, A, B, S, n, d, p_cov, tau)
    X = reshape(x(1:n*d), n, d);
    Z = reshape(x(n*d+1:end), p_cov, d);

    Y = X * S;
    S_mat = X * (S * X');
    S_mat(1:n+1:end) = 0;

    [psi_val, Psi_val, dpsi_val] = psi_functions(S_mat, tau);

    net_obj = sum((A(:) - S_mat(:)) .* psi_val(:) + Psi_val(:));

    B_pred = Z * X';
    cov_obj = -0.5 * sum((B(:) - B_pred(:)).^2);

    f = -(net_obj + cov_obj);

    if nargout > 1
        W_net = (A - S_mat) .* dpsi_val;
        W_net(1:n+1:end) = 0;

        grad_X_net = -2 * W_net * Y * S;

        resid_cov = B - B_pred;
        grad_X_cov = resid_cov' * Z;

        grad_X = -(grad_X_net + grad_X_cov);

        grad_Z = resid_cov * X;
        grad_Z = -grad_Z;

        g = [grad_X(:); grad_Z(:)];
    end
end
