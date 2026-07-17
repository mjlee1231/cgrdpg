% Compare trust-region vs quasi-newton algorithms for fminunc
% Test on single replication to see which finds better local minima

clear; clc;
addpath('core');

fprintf('Comparing fminunc Algorithms: trust-region vs quasi-newton\n');
fprintf('===========================================================\n\n');

% Parameters (match coverage test)
n = 1000;
p_cov = 500;
d = 3;
p = 2;
tau = 0.001;

% Set seed (rep 1)
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

% Generate network and covariates
A = double(rand(n) < P);
A = triu(A, 1);
A = A + A';
A(1:n+1:end) = 0;
B = Z0 * X0' + randn(p_cov, n);

fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P(:)), max(P(:)));

%% Test 1: Trust-Region (current approach)
fprintf('Test 1: Trust-Region Algorithm\n');
fprintf('-------------------------------\n');

options_tr = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');

t0 = tic;
[X_tr, ~, fval_tr, exitflag_tr, output_tr, ~] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options_tr);
time_tr = toc(t0);

[X_tr_aligned, ~] = procrustes_align(X_tr, X0);
sse_tr = sum((X_tr_aligned - X0).^2, 'all');

fprintf('Results:\n');
fprintf('  SSE: %.4f\n', sse_tr);
fprintf('  Time: %.1f sec\n', time_tr);
fprintf('  Exit flag: %d\n', exitflag_tr);
fprintf('  Outer iterations: %d\n\n', output_tr.iterations);

%% Test 2: Quasi-Newton (BFGS)
fprintf('Test 2: Quasi-Newton (BFGS) Algorithm\n');
fprintf('--------------------------------------\n');

options_qn = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true);

t0 = tic;
[X_qn, ~, fval_qn, exitflag_qn, output_qn, ~] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options_qn);
time_qn = toc(t0);

[X_qn_aligned, ~] = procrustes_align(X_qn, X0);
sse_qn = sum((X_qn_aligned - X0).^2, 'all');

fprintf('Results:\n');
fprintf('  SSE: %.4f\n', sse_qn);
fprintf('  Time: %.1f sec\n', time_qn);
fprintf('  Exit flag: %d\n', exitflag_qn);
fprintf('  Outer iterations: %d\n\n', output_qn.iterations);

%% Comparison
fprintf('Comparison\n');
fprintf('==========\n');
fprintf('SSE:  Trust-Region = %.4f,  Quasi-Newton = %.4f', sse_tr, sse_qn);
if sse_qn < sse_tr
    fprintf('  ← Quasi-Newton BETTER by %.1f%%\n', 100*(sse_tr - sse_qn)/sse_tr);
elseif sse_tr < sse_qn
    fprintf('  ← Trust-Region BETTER by %.1f%%\n', 100*(sse_qn - sse_tr)/sse_qn);
else
    fprintf('  (same)\n');
end

fprintf('Time: Trust-Region = %.1fs,  Quasi-Newton = %.1fs', time_tr, time_qn);
if time_qn < time_tr
    fprintf('  ← Quasi-Newton FASTER by %.1fx\n', time_tr/time_qn);
else
    fprintf('  ← Trust-Region FASTER by %.1fx\n', time_qn/time_tr);
end

fprintf('\nTarget (R rep 1): SSE = 8.23\n');
fprintf('Note: Even R''s Fisher scoring produces SSE 8-20, so anything under 30 is reasonable.\n');
