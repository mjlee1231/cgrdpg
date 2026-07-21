% Debug one-step cgrdpg estimator
% Check if update direction is correct

clear; clc;
addpath('core');

fprintf('Debugging One-Step cgrdpg Estimator\n');
fprintf('====================================\n\n');

% Full size test case (matching actual simulation)
n = 1000;
p_cov = 500;
d = 3;
p = 2;
tau = 0.001;
eps_clip = 1e-10;

% Set seed
rng(598 + 1);

% Generate data (OPTIMIZED design)
t = (1:n)' / n;
X0 = [0.42*t + 0.46, 0.27 * sin(2*pi*t) + 0.46, 0.20 * cos(4*pi*t)];
S = diag([1, 1, -1]);
Y0 = X0 * S;
Z0 = randn(p_cov, d);
P = X0 * Y0';

% Data
A = double(rand(n) < P);
A = triu(A, 1) + triu(A, 1)';
B = Z0 * X0' + randn(p_cov, n);

fprintf('True SSE baseline: 0.0000\n\n');

% ASE
[X_ase_unsigned, X_ase_signed, S_estimated] = fit_ase(A, d, p);
[X_ase_aligned, ~] = procrustes_align(X_ase_unsigned, X0);
sse_ase = sum((X_ase_aligned - X0).^2, 'all');
fprintf('ASE SSE: %.4f\n', sse_ase);

% Test SINGLE vertex update
i = 1;
x_i = X_ase_unsigned(i, :)';
XtX = X_ase_unsigned' * X_ase_unsigned;
Z_init = B * (X_ase_unsigned / XtX);

% Compute gradient and Fisher info
idx_j = setdiff(1:n, i);
Y_i = S_estimated * x_i;  % Y_i = S * x_i
s_i = X_ase_unsigned(idx_j, :) * Y_i;  % s_ij = x_j^T * S * x_i (linear predictor)
resid = A(i, idx_j)' - s_i;  % Residual (psi link treats s as pseudo-probability)
dpsi_val = dpsi(s_i, tau);  % Fisher weights using psi link

Y_j = X_ase_unsigned(idx_j, :) * S_estimated;  % FIXED: Y_j = X_j * S
grad_net = Y_j' * (resid .* dpsi_val);
grad_cov = Z_init' * (B(:, i) - Z_init * x_i);
grad_net_unscaled = grad_net;
grad_cov_unscaled = grad_cov;
grad = (grad_net + grad_cov) / (n + p_cov);

G_net = Y_j' * (Y_j .* dpsi_val);
G_cov = Z_init' * Z_init;
G = (G_net + G_cov) / (n + p_cov);

% Compute update direction
update_direction = (G + 1e-9 * eye(d)) \ grad;

fprintf('\nVertex %d analysis:\n', i);
fprintf('  Linear predictor s_ij: [%.6f, %.6f]\n', min(s_i), max(s_i));
fprintf('  dpsi weights range:    [%.2e, %.2e]\n', min(dpsi_val), max(dpsi_val));
fprintf('  Gradient_net (unscaled): %.6e\n', norm(grad_net_unscaled));
fprintf('  Gradient_cov (unscaled): %.6e\n', norm(grad_cov_unscaled));
fprintf('  Gradient (SCALED):       %.6e\n', norm(grad));
fprintf('  Update direction norm:   %.6f\n', norm(update_direction));
fprintf('  G condition number:      %.2e\n', cond(G));
fprintf('  Scaling factor (n+p):    %d\n', n + p_cov);

% Test both directions
x_plus = x_i + update_direction;
x_minus = x_i - update_direction;

X_temp_plus = X_ase_unsigned;
X_temp_plus(i, :) = x_plus';
[X_temp_plus_aligned, ~] = procrustes_align(X_temp_plus, X0);
dist_plus = sum((X_temp_plus_aligned - X0).^2, 'all');

X_temp_minus = X_ase_unsigned;
X_temp_minus(i, :) = x_minus';
[X_temp_minus_aligned, ~] = procrustes_align(X_temp_minus, X0);
dist_minus = sum((X_temp_minus_aligned - X0).^2, 'all');

fprintf('\nUpdate direction test:\n');
fprintf('  x_i + update:  distance = %.6f', dist_plus);
if dist_plus < sse_ase
    fprintf(' ✓ BETTER\n');
else
    fprintf(' ✗ WORSE\n');
end
fprintf('  x_i - update:  distance = %.6f', dist_minus);
if dist_minus < sse_ase
    fprintf(' ✓ BETTER\n');
else
    fprintf(' ✗ WORSE\n');
end

fprintf('\nConclusion:\n');
if dist_plus < dist_minus
    fprintf('  Should use: x_new = x + G^{-1} * grad (current)\n');
else
    fprintf('  Should use: x_new = x - G^{-1} * grad (OPPOSITE SIGN!)\n');
end

% Full one-step tests
fprintf('\n====================================\n');
fprintf('Full One-Step Test\n');
fprintf('====================================\n\n');

% Batch (fminunc) - true batch update
fprintf('Testing BATCH (fminunc)...\n');
[X_batch, ~] = compute_cgrdpg_one_step_batch_fminunc(A, B, X_ase_unsigned, S_estimated, tau);
[X_batch_aligned, ~] = procrustes_align(X_batch, X0);
sse_batch = sum((X_batch_aligned - X0).^2, 'all');

fprintf('BATCH (fminunc):\n');
fprintf('  SSE: %.4f', sse_batch);
if sse_batch < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n', 100*(sse_ase-sse_batch)/sse_ase);
else
    fprintf(' ✗ (%.1f%% worse)\n', 100*(sse_batch-sse_ase)/sse_ase);
end

% Jacobi - single Newton step
fprintf('\nTesting JACOBI (single Newton step)...\n');
[X_jacobi, ~] = compute_cgrdpg_one_step_jacobi(A, B, X_ase_unsigned, S_estimated, tau, 1.0);
[X_jacobi_aligned, ~] = procrustes_align(X_jacobi, X0);
sse_jacobi = sum((X_jacobi_aligned - X0).^2, 'all');

fprintf('JACOBI (single step):\n');
fprintf('  SSE: %.4f', sse_jacobi);
if sse_jacobi < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n\n', 100*(sse_ase-sse_jacobi)/sse_ase);
else
    fprintf(' ✗ (%.1f%% worse)\n\n', 100*(sse_jacobi-sse_ase)/sse_ase);
end

fprintf('\n====================================\n');
fprintf('COMPARISON\n');
fprintf('====================================\n');
fprintf('ASE SSE:             %.4f (baseline)\n', sse_ase);
fprintf('BATCH (fminunc):     %.4f', sse_batch);
if sse_batch < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n', 100*(sse_ase-sse_batch)/sse_ase);
else
    fprintf(' (%.1f%% worse)\n', 100*(sse_batch-sse_ase)/sse_ase);
end
fprintf('JACOBI (1 step):     %.4f', sse_jacobi);
if sse_jacobi < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n', 100*(sse_ase-sse_jacobi)/sse_ase);
else
    fprintf(' (%.1f%% worse)\n', 100*(sse_jacobi-sse_ase)/sse_ase);
end

fprintf('\nKey differences:\n');
fprintf('  BATCH uses fminunc: multiple inner iterations, all X updated simultaneously\n');
fprintf('  JACOBI uses 1 Newton step: updates computed using X_init, applied at once\n');
