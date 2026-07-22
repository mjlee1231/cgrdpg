function coverage_2d_no_info_single_rep_n500(rep_id)
% Vertex-wise Coverage: cgrdpg vs ASE vs OSE with NO COVARIATE INFORMATION
% Scenario: Z0 = 0, B contains only noise (no signal)
% This tests cgrdpg performance when covariates are uninformative
% n=500, p_cov=250, d=2 (RDPG)

if nargin < 1
    error('Usage: coverage_2d_no_info_single_rep_n500(rep_id)');
end

%% Setup paths
core_path = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), ...
    '3d_simulations', 'matlab', 'core');
addpath(core_path);

%% Parameters
n = 500;
p_cov = 250;
d = 2;
p = 2;  % RDPG: all positive eigenvalues
q = 0;
tau = 0.005;
eps_clip = 1e-10;
chi2_crit = chi2inv(0.95, d);

output_dir = 'results_2d_no_info_n500';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('============================================================================\n');
fprintf('  Vertex-wise Coverage: NO COVARIATE INFORMATION (n=%d)\n', n);
fprintf('  Replication %d/100\n', rep_id);
fprintf('  Z0 = 0, B = pure noise\n');
fprintf('  Methods: cgrdpg, ASE, OSE\n');
fprintf('  S = diag([1, 1]) (RDPG)\n');
fprintf('============================================================================\n\n');

rep_start = tic;

% Set seed for reproducibility
rng(598 + rep_id);

%% 1. Generate true latent positions (circular design)
theta = pi * (1:n)' / (n - 1);
X0 = zeros(n, d);
X0(:, 1) = 0.28 * sin(theta) + 0.42;
X0(:, 2) = 0.28 * cos(theta) + 0.42;

S = eye(d);  % RDPG: positive definite
Y0 = X0 * S;

% 2. NO COVARIATE INFORMATION: Z0 = 0
Z0 = zeros(p_cov, d);  % Zero signal matrix

% Edge probabilities
P = X0 * Y0';
fprintf('Edge probability range: [%.4f, %.4f]\n\n', min(P(:)), max(P(:)));

%% 2. Generate data
A = double(rand(n) < P);
A = triu(A, 1);
A = A + A';  % Symmetric
A(1:n+1:end) = 0;  % No self-loops

% B contains ONLY NOISE (no signal from Z0 * X0')
B = randn(p_cov, n);

%% 3. ASE (computed first to initialize cgrdpg)
fprintf('Computing ASE...\n');
t0 = tic;
[X_ase_unsigned, X_ase_signed, S_estimated, ase_eigenvalues, ase_eigenvectors] = fit_ase(A, d, p);
ase_time = toc(t0);
[X_ase, ~] = procrustes_align(X_ase_unsigned, X0);
fprintf('ASE: time=%.1fs, S_est=diag([%+d,%+d])\n', ...
    ase_time, diag(S_estimated));

%% 4. Fit cgrdpg using fminunc (for better convergence)
fprintf('Fitting cgrdpg with fminunc...\n');
t0 = tic;

% Use fminunc with surrogate objective
options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');

% Initialize from ASE
X_init = X_ase_unsigned;
Y_init = X_init * S_estimated;
Z_init = B * X_init / (X_init' * X_init);

x0 = [X_init(:); Y_init(:); Z_init(:)];

% Objective function
obj_fun = @(x) surrogate_objective_2d(x, A, B, n, p_cov, d, tau);

% Optimize
[x_opt, fval, exitflag, output] = fminunc(obj_fun, x0, options);

cgrdpg_time = toc(t0);

% Extract solution
X_cgrdpg = reshape(x_opt(1:n*d), n, d);
Y_cgrdpg = reshape(x_opt(n*d+1:2*n*d), n, d);
Z_cgrdpg = reshape(x_opt(2*n*d+1:end), p_cov, d);

% Align to true positions
[X_cgrdpg, ~] = procrustes_align(X_cgrdpg, X0);
Y_cgrdpg = X_cgrdpg * S_estimated;
Z_cgrdpg = B * X_cgrdpg / (X_cgrdpg' * X_cgrdpg);

fprintf('cgrdpg: exitflag=%d, iters=%d, time=%.1fs, fval=%.4f\n', ...
    exitflag, output.iterations, cgrdpg_time, fval);

%% 5. OSE (one-step from ASE)
fprintf('Computing OSE...\n');
t0 = tic;
X_ose_unsigned = compute_ose_step(A, X_ase_unsigned, eps_clip);
ose_step_time = toc(t0);
ose_time = ase_time + ose_step_time;
[X_ose, ~] = procrustes_align(X_ose_unsigned, X0);
fprintf('OSE: time=%.1fs (ASE: %.1fs + step: %.1fs)\n', ...
    ose_time, ase_time, ose_step_time);

%% 6. SSE
sse_cgrdpg = sum((X_cgrdpg - X0).^2, 'all');
sse_ase = sum((X_ase - X0).^2, 'all');
sse_ose = sum((X_ose - X0).^2, 'all');

fprintf('SSE: cgrdpg=%.4f  ASE=%.4f  OSE=%.4f\n\n', sse_cgrdpg, sse_ase, sse_ose);

%% 7. Vertex-wise coverage
fprintf('Computing vertex-wise coverage for all %d vertices...\n', n);

results = struct();
results.cgrdpg_true = nan(n, 1);
results.cgrdpg_plugin = nan(n, 1);
results.ase_true = nan(n, 1);
results.ase_plugin = nan(n, 1);
results.ose_true = nan(n, 1);
results.ose_plugin = nan(n, 1);

t0 = tic;
for i = 1:n
    if mod(i, 100) == 0
        fprintf('  Vertex %d/%d\n', i, n);
    end

    % cgrdpg
    results.cgrdpg_true(i) = check_coverage(...
        X0(i,:) - X_cgrdpg(i,:), ...
        compute_fisher_info_cgrdpg(i, X0, Y0, Z0, tau), ...
        n + p_cov, chi2_crit);

    results.cgrdpg_plugin(i) = check_coverage(...
        X0(i,:) - X_cgrdpg(i,:), ...
        compute_fisher_info_cgrdpg(i, X_cgrdpg, Y_cgrdpg, Z_cgrdpg, tau), ...
        n + p_cov, chi2_crit);

    % ASE
    results.ase_true(i) = check_coverage(...
        X_ase(i,:) - X0(i,:), ...
        compute_precision_ase(i, X0, S, eps_clip), ...
        1, chi2_crit);

    results.ase_plugin(i) = check_coverage(...
        X_ase(i,:) - X0(i,:), ...
        compute_precision_ase(i, X_ase, S_estimated, eps_clip), ...
        1, chi2_crit);

    % OSE
    results.ose_true(i) = check_coverage(...
        X_ose(i,:) - X0(i,:), ...
        compute_precision_ose(i, X0, eps_clip), ...
        1, chi2_crit);

    results.ose_plugin(i) = check_coverage(...
        X_ose(i,:) - X0(i,:), ...
        compute_precision_ose(i, X_ose, eps_clip), ...
        1, chi2_crit);
end

cov_time = toc(t0);

% Overall coverage
overall_cov = struct();
overall_cov.cgrdpg_true = mean(results.cgrdpg_true, 'omitnan');
overall_cov.cgrdpg_plugin = mean(results.cgrdpg_plugin, 'omitnan');
overall_cov.ase_true = mean(results.ase_true, 'omitnan');
overall_cov.ase_plugin = mean(results.ase_plugin, 'omitnan');
overall_cov.ose_true = mean(results.ose_true, 'omitnan');
overall_cov.ose_plugin = mean(results.ose_plugin, 'omitnan');

% Count NAs
n_na = struct();
n_na.cgrdpg_true = sum(isnan(results.cgrdpg_true));
n_na.cgrdpg_plugin = sum(isnan(results.cgrdpg_plugin));
n_na.ase_true = sum(isnan(results.ase_true));
n_na.ase_plugin = sum(isnan(results.ase_plugin));
n_na.ose_true = sum(isnan(results.ose_true));
n_na.ose_plugin = sum(isnan(results.ose_plugin));

rep_time = toc(rep_start) / 60;  % minutes

fprintf('\nOverall coverage (this rep):\n');
fprintf('  cgrdpg-TRUE      %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.cgrdpg_true, n_na.cgrdpg_true);
fprintf('  cgrdpg-PLUGIN    %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.cgrdpg_plugin, n_na.cgrdpg_plugin);
fprintf('  ASE-TRUE         %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.ase_true, n_na.ase_true);
fprintf('  ASE-PLUGIN       %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.ase_plugin, n_na.ase_plugin);
fprintf('  OSE-TRUE         %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.ose_true, n_na.ose_true);
fprintf('  OSE-PLUGIN       %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.ose_plugin, n_na.ose_plugin);
fprintf('\nTotal rep time: %.2f min\n', rep_time);

%% 8. Save results
out_file = fullfile(output_dir, sprintf('rep_%03d.mat', rep_id));
save(out_file, ...
    'rep_id', 'n', 'p_cov', 'd', 'tau', 'S', 'S_estimated', ...
    'results', 'overall_cov', 'n_na', ...
    'X0', 'X_cgrdpg', 'Y_cgrdpg', 'Z_cgrdpg', 'X_ase', 'X_ose', ...
    'sse_cgrdpg', 'sse_ase', 'sse_ose', ...
    'cgrdpg_time', 'ase_time', 'ose_time', 'cov_time', 'rep_time', ...
    'exitflag', 'output', ...
    'ase_eigenvalues', 'ase_eigenvectors');

fprintf('\nResults saved to: %s\n', out_file);
fprintf('============================================================================\n');

end
