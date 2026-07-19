% Test new latent position design vs current on single replication
% Compare optimization performance with better eigenvalue separation

clear; clc;
addpath('core');

fprintf('Comparing Old vs New Latent Position Design\n');
fprintf('============================================\n\n');

% Parameters
n = 1000;
p_cov = 500;
d = 3;
p = 2;
tau = 0.001;

% Set seed
rng(598 + 1);

t = (1:n)' / n;
S = diag([1, 1, -1]);

%% Design 1: Current (degenerate eigenvalues)
fprintf('Design 1: Current 3D Helix\n');
fprintf('---------------------------\n');

X0_old = [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)];

Y0_old = X0_old * S;
Z0_old = randn(p_cov, d);
P_old = X0_old * Y0_old';

fprintf('Edge prob range: [%.4f, %.4f]\n', min(P_old(:)), max(P_old(:)));

eigvals_old = eig(P_old);
[~, idx] = sort(abs(eigvals_old), 'descend');
top3_old = eigvals_old(idx(1:3));
fprintf('Top 3 eigenvalues: %+.2f, %+.2f, %+.2f\n', top3_old);
fprintf('Ratios: λ₂/λ₁=%.3f, λ₃/λ₂=%.3f\n', ...
    abs(top3_old(2))/abs(top3_old(1)), abs(top3_old(3))/abs(top3_old(2)));
fprintf('Gap 2-3: %.4f ← PROBLEM: Nearly zero!\n\n', ...
    abs(abs(top3_old(2)) - abs(top3_old(3))));

%% Design 2: New (well-separated eigenvalues)
fprintf('Design 2: Linear + Periodic\n');
fprintf('----------------------------\n');

X0_new = [0.3*t + 0.5, ...
          0.15 * sin(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)];

Y0_new = X0_new * S;
Z0_new = Z0_old;  % Same Z for fair comparison
P_new = X0_new * Y0_new';

fprintf('Edge prob range: [%.4f, %.4f]\n', min(P_new(:)), max(P_new(:)));

eigvals_new = eig(P_new);
[~, idx] = sort(abs(eigvals_new), 'descend');
top3_new = eigvals_new(idx(1:3));
fprintf('Top 3 eigenvalues: %+.2f, %+.2f, %+.2f\n', top3_new);
fprintf('Ratios: λ₂/λ₁=%.3f, λ₃/λ₂=%.3f\n', ...
    abs(top3_new(2))/abs(top3_new(1)), abs(top3_new(3))/abs(top3_new(2)));
fprintf('Gap 2-3: %.4f ← BETTER: Well separated!\n\n', ...
    abs(abs(top3_new(2)) - abs(top3_new(3))));

%% Generate data for both
% Old design
A_old = double(rand(n) < P_old);
A_old = triu(A_old, 1);
A_old = A_old + A_old';
A_old(1:n+1:end) = 0;
B_old = Z0_old * X0_old' + randn(p_cov, n);

% New design
rng(598 + 1);  % Reset to get same random graph structure
A_new = double(rand(n) < P_new);
A_new = triu(A_new, 1);
A_new = A_new + A_new';
A_new(1:n+1:end) = 0;
B_new = Z0_new * X0_new' + randn(p_cov, n);

%% Test batch surrogate on both
fprintf('========================================\n');
fprintf('Test: Batch Surrogate Optimization\n');
fprintf('========================================\n\n');

options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');

% Old design
fprintf('Old design (degenerate eigenvalues):\n');
t0 = tic;
[X_opt_old, ~, ~, ~, output_old, ~] = ...
    fit_grdpg_fminunc_surrogate(A_old, B_old, d, p, tau, options);
time_old = toc(t0);

[X_aligned_old, ~] = procrustes_align(X_opt_old, X0_old);
sse_old = sum((X_aligned_old - X0_old).^2, 'all');

fprintf('  SSE: %.4f\n', sse_old);
fprintf('  Time: %.1f sec\n', time_old);
fprintf('  Iterations: %d\n\n', output_old.iterations);

% New design
fprintf('New design (well-separated eigenvalues):\n');
t0 = tic;
[X_opt_new, ~, ~, ~, output_new, ~] = ...
    fit_grdpg_fminunc_surrogate(A_new, B_new, d, p, tau, options);
time_new = toc(t0);

[X_aligned_new, ~] = procrustes_align(X_opt_new, X0_new);
sse_new = sum((X_aligned_new - X0_new).^2, 'all');

fprintf('  SSE: %.4f\n', sse_new);
fprintf('  Time: %.1f sec\n', time_new);
fprintf('  Iterations: %d\n\n', output_new.iterations);

%% Compare
fprintf('========================================\n');
fprintf('Comparison\n');
fprintf('========================================\n\n');

fprintf('SSE:\n');
fprintf('  Old design: %.4f\n', sse_old);
fprintf('  New design: %.4f', sse_new);
if sse_new < sse_old
    improvement = 100 * (sse_old - sse_new) / sse_old;
    fprintf(' ✓ BETTER by %.1f%%\n', improvement);
else
    degradation = 100 * (sse_new - sse_old) / sse_old;
    fprintf(' ✗ Worse by %.1f%%\n', degradation);
end

fprintf('\nConvergence:\n');
fprintf('  Old design: %d iterations\n', output_old.iterations);
fprintf('  New design: %d iterations\n', output_new.iterations);

fprintf('\n========================================\n');
if sse_new < sse_old
    fprintf('CONCLUSION: New design improves optimization!\n');
    fprintf('Recommend: Switch to new latent position design\n');
else
    fprintf('CONCLUSION: New design does not improve optimization\n');
    fprintf('Recommend: Keep current design or investigate further\n');
end
fprintf('========================================\n');
