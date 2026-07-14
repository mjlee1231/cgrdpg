% Test vectorization and reshaping correctness
clear; clc;

% Add core folder to path
addpath('core');

fprintf('Testing Vectorization and Index Ordering\n');
fprintf('========================================\n\n');

%% Very small problem to manually verify
n = 3;
p_cov = 2;
d = 2;
tau = 0.001;
rng(123);

% Generate tiny data
theta = pi * (1:n)' / (n - 1);
X_true = [0.28 * sin(theta) + 0.42, ...
          0.28 * cos(theta) + 0.42];
Z_true = randn(p_cov, d) * 0.1;
S = eye(d);
Y_true = X_true * S;
P_net = X_true * Y_true';
A = double(rand(n) < P_net);
A = triu(A, 1);
A = A + A';
A(1:n+1:end) = 0;
B = Z_true * X_true' + randn(p_cov, n) * 0.1;

% Initialize
X_init = X_true + randn(n, d) * 0.01;
S_estimated = eye(d);
Y_hat = X_init * S_estimated;
Z_hat = (X_init \ B')';

fprintf('Problem size: n=%d, d=%d, p_cov=%d\n\n', n, d, p_cov);
fprintf('X_init (n x d):\n');
disp(X_init);

% Pack into vector (column-major)
x0 = X_init(:);
fprintf('x0 = X_init(:) vectorization:\n');
fprintf('  x0(1) = X_init(1,1) = %.6f\n', x0(1));
fprintf('  x0(2) = X_init(2,1) = %.6f\n', x0(2));
fprintf('  x0(3) = X_init(3,1) = %.6f\n', x0(3));
fprintf('  x0(4) = X_init(1,2) = %.6f\n', x0(4));
fprintf('  x0(5) = X_init(2,2) = %.6f\n', x0(5));
fprintf('  x0(6) = X_init(3,2) = %.6f\n\n', x0(6));

% Define objective
objective = @(x) surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau);

% Compute analytical gradient
[f0, g_analytical] = objective(x0);

fprintf('Analytical gradient g_analytical:\n');
fprintf('  g(1) [node 1, dim 1]: %.6e\n', g_analytical(1));
fprintf('  g(2) [node 2, dim 1]: %.6e\n', g_analytical(2));
fprintf('  g(3) [node 3, dim 1]: %.6e\n', g_analytical(3));
fprintf('  g(4) [node 1, dim 2]: %.6e\n', g_analytical(4));
fprintf('  g(5) [node 2, dim 2]: %.6e\n', g_analytical(5));
fprintf('  g(6) [node 3, dim 2]: %.6e\n\n', g_analytical(6));

%% Manual finite difference for specific components
eps_fd = 1e-7;

fprintf('Manual finite difference check:\n');
fprintf('========================================\n\n');

% Check x(1) = X(1,1) [node 1, dimension 1]
fprintf('Component 1: X(1,1) [node 1, dim 1]\n');
x_plus = x0; x_plus(1) = x_plus(1) + eps_fd;
x_minus = x0; x_minus(1) = x_minus(1) - eps_fd;
fd_1 = (objective(x_plus) - objective(x_minus)) / (2 * eps_fd);
fprintf('  Analytical: %.6e\n', g_analytical(1));
fprintf('  Finite diff: %.6e\n', fd_1);
fprintf('  Ratio: %.6f\n\n', g_analytical(1) / fd_1);

% Check x(4) = X(1,2) [node 1, dimension 2]
fprintf('Component 4: X(1,2) [node 1, dim 2]\n');
x_plus = x0; x_plus(4) = x_plus(4) + eps_fd;
x_minus = x0; x_minus(4) = x_minus(4) - eps_fd;
fd_4 = (objective(x_plus) - objective(x_minus)) / (2 * eps_fd);
fprintf('  Analytical: %.6e\n', g_analytical(4));
fprintf('  Finite diff: %.6e\n', fd_4);
fprintf('  Ratio: %.6f\n\n', g_analytical(4) / fd_4);

% Check x(2) = X(2,1) [node 2, dimension 1]
fprintf('Component 2: X(2,1) [node 2, dim 1]\n');
x_plus = x0; x_plus(2) = x_plus(2) + eps_fd;
x_minus = x0; x_minus(2) = x_minus(2) - eps_fd;
fd_2 = (objective(x_plus) - objective(x_minus)) / (2 * eps_fd);
fprintf('  Analytical: %.6e\n', g_analytical(2));
fprintf('  Finite diff: %.6e\n', fd_2);
fprintf('  Ratio: %.6f\n\n', g_analytical(2) / fd_2);

%% Check reshape consistency
fprintf('========================================\n');
fprintf('Reshaping consistency check:\n');
fprintf('========================================\n\n');

X_reshaped = reshape(x0, n, d);
fprintf('reshape(x0, n, d):\n');
disp(X_reshaped);

fprintf('Does reshape(x0, n, d) == X_init? %d\n', isequal(X_reshaped, X_init));
fprintf('Max difference: %.6e\n\n', max(abs(X_reshaped(:) - X_init(:))));

% Reverse check
grad_X_manual = reshape(g_analytical, n, d);
fprintf('Gradient as matrix (reshape(g, n, d)):\n');
disp(grad_X_manual);

%% Nested function
function [f, g] = surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau)
    % Same as in test_gradient.m
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
