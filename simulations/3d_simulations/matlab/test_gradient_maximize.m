% Test gradient correctness - MAXIMIZING version
clear; clc;

% Add core folder to path
addpath('core');

fprintf('Testing Gradient Correctness (MAXIMIZING)\n');
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
eigvals_sorted = eigvals(idx);
X_init = V(:, idx(1:d)) * diag(sqrt(abs(eigvals_sorted(1:d))));

% Construct signature from p parameter (p positive, q=d-p negative)
q = d - p;
S_estimated = diag([ones(p, 1); -ones(q, 1)]);

fprintf('Signature: S = diag([');
fprintf('%+d ', diag(S_estimated)');
fprintf('])\n\n');

% Fix Y_hat and Z_hat (as in surrogate algorithm)
Y_hat = X_init * S_estimated;
% Stably solve Z_hat from: Z_hat * X_init' = B
Z_hat = (X_init \ B')';

% Pack X only (Y_hat and Z_hat are FIXED)
x0 = X_init(:);

% Define objective with FIXED Y_hat and Z_hat (MAXIMIZING VERSION)
objective = @(x) surrogate_objective_gradient_maximize(x, A, B, Y_hat, Z_hat, n, d, tau);

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

%% MAXIMIZING version - no negation
function [f, g] = surrogate_objective_gradient_maximize(x, A, B, Y_hat, Z_hat, n, d, tau)
    % Compute surrogate objective and gradient (MAXIMIZING version)
    % Using logical masking to cleanly exclude self-loops

    % 1. Unpack X (column-major)
    X = reshape(x, n, d);

    % 2. Network probabilities (keep diagonal values, don't set to 0)
    S_mat = X * Y_hat';

    % Create off-diagonal mask
    is_off_diag = ~eye(n);

    % 3. Compute psi and Psi values
    [psi_val, Psi_val, dpsi_val] = psi_functions(S_mat, tau);

    % 4. Network component calculation
    % Use masking to exclude diagonal (self-loop) from summation
    net_obj = sum(A(is_off_diag) .* psi_val(is_off_diag) + Psi_val(is_off_diag));

    % 5. Covariate component: -0.5 * ||B - Z_hat*X'||_F^2
    B_pred = Z_hat * X';
    cov_obj = -0.5 * sum((B(:) - B_pred(:)).^2);

    % 6. Total objective (MAXIMIZING - no negation)
    f = (net_obj + cov_obj);

    % 7. Compute gradient if requested
    if nargout > 1
        % (1) Compute full weight matrix (including diagonal)
        W_full = (A - S_mat) .* dpsi_val;

        % (2) Network gradient: Y_hat is fixed constant, coefficient is 1
        % MAXIMIZING - no negation
        grad_X_net = W_full * Y_hat;

        % (3) Correction: Remove diagonal contribution that was incorrectly included
        diag_W = diag(W_full);
        grad_X_net = grad_X_net - diag_W .* Y_hat;

        % (4) Covariate gradient (MAXIMIZING - no negation)
        resid_cov = B - B_pred;
        grad_X_cov = resid_cov' * Z_hat;

        % (5) Total gradient
        grad_X = grad_X_net + grad_X_cov;

        % Pack gradient (column-major)
        g = grad_X(:);
    end
end
