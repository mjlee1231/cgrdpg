% Debug one-step cgrdpg estimator
% Check if update direction is correct

clear; clc;
addpath('core');

fprintf('Debugging One-Step cgrdpg Estimator\n');
fprintf('====================================\n\n');

% Small test case
n = 100;
p_cov = 50;
d = 3;
p = 2;
tau = 0.001;
eps_clip = 1e-10;

% Set seed
rng(598 + 1);

% Generate data (NEW design)
t = (1:n)' / n;
X0 = [0.3*t + 0.5, 0.15 * sin(2*pi*t) + 0.6, 0.1 * cos(4*pi*t)];
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
Y_init = X_ase_unsigned * S_estimated;
XtX = X_ase_unsigned' * X_ase_unsigned;
Z_init = B * (X_ase_unsigned / XtX);

% Compute gradient and Fisher info
idx_j = setdiff(1:n, i);
s_i = Y_init(idx_j, :) * x_i;
p_i = max(min(1 ./ (1 + exp(-s_i / tau)), 1 - eps_clip), eps_clip);
resid = A(i, idx_j)' - p_i;
dpsi_val = 1 ./ (tau * p_i .* (1 - p_i));

grad_net = Y_init(idx_j, :)' * (resid .* dpsi_val);
grad_cov = Z_init' * (B(:, i) - Z_init * x_i);
grad = (grad_net + grad_cov) / (n + p_cov);

G_net = Y_init(idx_j, :)' * (Y_init(idx_j, :) .* dpsi_val);
G_cov = Z_init' * Z_init;
G = (G_net + G_cov) / (n + p_cov);

% Compute update direction
update_direction = (G + 1e-9 * eye(d)) \ grad;

fprintf('\nVertex %d analysis:\n', i);
fprintf('  Gradient norm:        %.6f\n', norm(grad));
fprintf('  Update direction norm: %.6f\n', norm(update_direction));
fprintf('  G condition number:    %.2e\n', cond(G));

% Test both directions
x_plus = x_i + update_direction;
x_minus = x_i - update_direction;

[~, ~, dist_plus] = procrustes_align([x_plus'; X_ase_unsigned(2:end,:)], X0);
[~, ~, dist_minus] = procrustes_align([x_minus'; X_ase_unsigned(2:end,:)], X0);

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

% Full one-step with BOTH signs
fprintf('\n====================================\n');
fprintf('Full One-Step Test\n');
fprintf('====================================\n\n');

% Current implementation (+)
[X_plus, ~] = compute_cgrdpg_one_step_batch(A, B, X_ase_unsigned, S_estimated, tau);
[X_plus_aligned, ~] = procrustes_align(X_plus, X0);
sse_plus = sum((X_plus_aligned - X0).^2, 'all');

fprintf('Current (x + update):\n');
fprintf('  SSE: %.4f', sse_plus);
if sse_plus < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n', 100*(sse_ase-sse_plus)/sse_ase);
else
    fprintf(' ✗ (%.1f%% worse)\n\n', 100*(sse_plus-sse_ase)/sse_ase);
end

% Try opposite sign (-)
X_minus = zeros(n, d);
for i = 1:n
    x_i = X_ase_unsigned(i, :)';
    idx_j = setdiff(1:n, i);
    s_i = Y_init(idx_j, :) * x_i;
    p_i = max(min(1 ./ (1 + exp(-s_i / tau)), 1 - eps_clip), eps_clip);
    resid = A(i, idx_j)' - p_i;
    dpsi_val = 1 ./ (tau * p_i .* (1 - p_i));

    grad_net = Y_init(idx_j, :)' * (resid .* dpsi_val);
    grad_cov = Z_init' * (B(:, i) - Z_init * x_i);
    grad = (grad_net + grad_cov) / (n + p_cov);

    G_net = Y_init(idx_j, :)' * (Y_init(idx_j, :) .* dpsi_val);
    G_cov = Z_init' * Z_init;
    G = (G_net + G_cov) / (n + p_cov);

    % OPPOSITE SIGN
    X_minus(i, :) = (x_i - (G + 1e-9 * eye(d)) \ grad)';
end

[X_minus_aligned, ~] = procrustes_align(X_minus, X0);
sse_minus = sum((X_minus_aligned - X0).^2, 'all');

fprintf('Opposite (x - update):\n');
fprintf('  SSE: %.4f', sse_minus);
if sse_minus < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n', 100*(sse_ase-sse_minus)/sse_ase);
else
    fprintf(' ✗ (%.1f%% worse)\n\n', 100*(sse_minus-sse_ase)/sse_ase);
end

fprintf('\n====================================\n');
fprintf('RECOMMENDATION\n');
fprintf('====================================\n');
if sse_minus < sse_plus && sse_minus < sse_ase
    fprintf('✓ Use OPPOSITE SIGN: x_new = x - G^{-1} * grad\n');
    fprintf('  This will improve SSE from %.1f to %.1f\n', sse_ase, sse_minus);
elseif sse_plus < sse_ase
    fprintf('✓ Current sign is correct: x_new = x + G^{-1} * grad\n');
    fprintf('  SSE improves from %.1f to %.1f\n', sse_ase, sse_plus);
else
    fprintf('✗ One-step does not improve over ASE!\n');
    fprintf('  Check gradient/Fisher info calculation\n');
end
