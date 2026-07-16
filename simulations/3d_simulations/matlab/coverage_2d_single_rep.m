function coverage_2d_single_rep(rep_id)
% COVERAGE_2D_SINGLE_REP Compute vertex-wise coverage for 2D GRDPG
%
% Single replication for cgrdpg vertex-wise coverage testing
% Tests: cgrdpg-TRUE and cgrdpg-PLUGIN methods
%
% Inputs:
%   rep_id - replication number (1-100)
%
% Parameters:
%   n = 500, p_cov = 250, d = 2
%   S = eye(2) (positive definite)
%   Latent curve: semicircle in 2D

% Add core folder to path
addpath('core');

% Parameters
n = 500;
p_cov = 250;
d = 2;
p = 2;  % All positive (identity signature)
maxit = 30;
tol = 0.01;
tau = 0.001;
eps_clip = 1e-10;
chi2_crit = chi2inv(0.95, d);

% Output directory
output_dir = 'results_matlab_2d_coverage';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('============================================================================\n');
fprintf('  cgrdpg Vertex-wise Coverage: 2D GRDPG, n=%d\n', n);
fprintf('  Replication %d/100\n', rep_id);
fprintf('  Methods: cgrdpg-TRUE, cgrdpg-PLUGIN\n');
fprintf('  S = eye(2), p_cov=%d\n', p_cov);
fprintf('============================================================================\n\n');

rep_start = tic;

% Set seed for reproducibility
rng(598 + rep_id);

%% 1. Generate true latent positions (2D semicircle)
theta = pi * (1:n)' / (n - 1);
X0 = [0.28 * sin(theta) + 0.42, ...
      0.28 * cos(theta) + 0.42];

S = eye(d);  % Identity signature (all positive)
Y0 = X0 * S;  % For positive definite, Y0 = X0
Z0 = randn(p_cov, d);

% Edge probabilities
P = X0 * Y0';
fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P(:)), max(P(:)));

%% 2. Generate data
A = double(rand(n) < P);
A = triu(A, 1);
A = A + A';  % Symmetric
A(1:n+1:end) = 0;  % No self-loops

B = Z0 * X0' + randn(p_cov, n);

%% 3. Fit cgrdpg
fprintf('Fitting cgrdpg...\n');
t0 = tic;

options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');

[X_opt, Z_opt, fval, exitflag, output] = ...
    fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options);

cgrdpg_time = toc(t0);

% Procrustes alignment
[X_cgrdpg, ~] = procrustes_align(X_opt, X0);
Y_cgrdpg = X_cgrdpg * S;
% Z_cgrdpg = B * X_cgrdpg * inv(X_cgrdpg' * X_cgrdpg)
Z_cgrdpg = (X_cgrdpg \ B')';

fprintf('cgrdpg: converged=%d, iters=%d, time=%.1fs\n', ...
    exitflag, output.iterations, cgrdpg_time);

% SSE
sse_cgrdpg = sum((X_cgrdpg - X0).^2, 'all');
fprintf('SSE cgrdpg=%.4f\n\n', sse_cgrdpg);

%% 4. Vertex-wise coverage
fprintf('Computing vertex-wise coverage for all %d vertices...\n', n);

results = struct();
results.cgrdpg_true = nan(n, 1);
results.cgrdpg_plugin = nan(n, 1);

t0 = tic;
for i = 1:n
    if mod(i, 100) == 0
        fprintf('  Vertex %d/%d\n', i, n);
    end

    % Error vector
    err = X0(i,:) - X_cgrdpg(i,:);

    % cgrdpg-TRUE: use true X0, Y0, Z0
    G_true = compute_fisher_info_cgrdpg(i, X0, Y0, Z0, tau);
    results.cgrdpg_true(i) = check_coverage(err, G_true, chi2_crit, n + p_cov);

    % cgrdpg-PLUGIN: use estimated X_cgrdpg, Y_cgrdpg, Z_cgrdpg
    G_plugin = compute_fisher_info_cgrdpg(i, X_cgrdpg, Y_cgrdpg, Z_cgrdpg, tau);
    results.cgrdpg_plugin(i) = check_coverage(err, G_plugin, chi2_crit, n + p_cov);
end

cov_time = toc(t0);

% Overall coverage rates
overall_cov = struct();
overall_cov.cgrdpg_true = mean(results.cgrdpg_true, 'omitnan');
overall_cov.cgrdpg_plugin = mean(results.cgrdpg_plugin, 'omitnan');

% Count NAs
n_na = struct();
n_na.cgrdpg_true = sum(isnan(results.cgrdpg_true));
n_na.cgrdpg_plugin = sum(isnan(results.cgrdpg_plugin));

rep_time = toc(rep_start) / 60;  % minutes

fprintf('\nOverall coverage (this rep):\n');
fprintf('  cgrdpg-TRUE      %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.cgrdpg_true, n_na.cgrdpg_true);
fprintf('  cgrdpg-PLUGIN    %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.cgrdpg_plugin, n_na.cgrdpg_plugin);
fprintf('\nTotal rep time: %.2f min\n', rep_time);

%% 5. Save results
out_file = fullfile(output_dir, sprintf('rep_%03d.mat', rep_id));
save(out_file, ...
    'rep_id', 'n', 'p_cov', 'd', 'tau', 'S', ...
    'results', 'overall_cov', 'n_na', 'sse_cgrdpg', ...
    'cgrdpg_time', 'cov_time', 'rep_time', ...
    'exitflag', 'output');

fprintf('\nResults saved to: %s\n', out_file);
fprintf('============================================================================\n');

end
