function coverage_3d_single_rep_jacobi(rep_id)
% COVERAGE_3D_SINGLE_REP_JACOBI Vertex-wise coverage using Jacobi method
%
% Same as coverage_3d_single_rep but uses Jacobi coordinate ascent
% instead of batch surrogate for cgrdpg optimization
%
% Inputs:
%   rep_id - replication number (1-100)

% Add core folder to path
addpath('core');

% Parameters
n = 1000;
p_cov = 500;
d = 3;
p = 2;  % 2 positive, 1 negative
q = 1;
maxit = 30;
tol = 0.005;  % Match batch surrogate
tau = 0.001;
eps_clip = 1e-10;
chi2_crit = chi2inv(0.95, d);

% Output directory
output_dir = 'results_matlab_3d_coverage_jacobi';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('============================================================================\n');
fprintf('  ASE/OSE/cgrdpg Vertex-wise Coverage: 3D GRDPG, n=%d (JACOBI METHOD)\n', n);
fprintf('  Replication %d/100\n', rep_id);
fprintf('  Methods: cgrdpg-TRUE/PLUGIN, ASE-TRUE/PLUGIN, OSE-TRUE/PLUGIN\n');
fprintf('  S = diag([1, 1, -1]), p_cov=%d\n', p_cov);
fprintf('============================================================================\n\n');

rep_start = tic;

% Set seed for reproducibility
rng(598 + rep_id);

%% 1. Generate true latent positions (linear + periodic, well-separated eigenvalues)
t = (1:n)' / n;
X0 = [0.3*t + 0.5, ...
      0.15 * sin(2*pi*t) + 0.6, ...
      0.1 * cos(4*pi*t)];

S = diag([1, 1, -1]);  % Indefinite signature
Y0 = X0 * S;
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

%% 3. ASE (computed first to get estimated signature)
fprintf('Computing ASE...\n');
t0 = tic;
[X_ase_unsigned, X_ase_signed, S_estimated] = fit_ase(A, d, p);
ase_time = toc(t0);
[X_ase, ~] = procrustes_align(X_ase_unsigned, X0);
fprintf('ASE: time=%.1fs, S_est=diag([%+d,%+d,%+d])\n', ...
    ase_time, diag(S_estimated));

%% 4. Fit cgrdpg using JACOBI coordinate ascent
fprintf('Fitting cgrdpg with Jacobi coordinate ascent...\n');
t0 = tic;

% Use Jacobi with adaptive step size (best from tuning: 0.01 initial)
step_size_init = 0.01;
[X_opt, Z_opt, fval, exitflag, output, ~] = ...
    fit_grdpg_jacobi_custom_step(A, B, d, p, tau, maxit, tol, step_size_init);

cgrdpg_time = toc(t0);

% Procrustes alignment
[X_cgrdpg, ~] = procrustes_align(X_opt, X0);
% CRITICAL: Use estimated signature S_estimated (not true S) for PLUGIN method
Y_cgrdpg = X_cgrdpg * S_estimated;
Z_cgrdpg = (X_cgrdpg \ B')';

fprintf('cgrdpg (Jacobi): converged=%d, iters=%d, time=%.1fs\n', ...
    exitflag, output.iterations, cgrdpg_time);

%% 5. OSE
fprintf('Computing OSE...\n');
t0 = tic;
X_ose_raw = compute_ose_step(A, X_ase_unsigned, X_ase_signed, eps_clip);
ose_step_time = toc(t0);
ose_time = ase_time + ose_step_time;
[X_ose, ~] = procrustes_align(X_ose_raw, X0);
fprintf('OSE: time=%.1fs (ASE:%.1fs + step:%.1fs)\n', ...
    ose_time, ase_time, ose_step_time);

% SSE
sse_cgrdpg = sum((X_cgrdpg - X0).^2, 'all');
sse_ase = sum((X_ase - X0).^2, 'all');
sse_ose = sum((X_ose - X0).^2, 'all');
fprintf('SSE: cgrdpg=%.4f  ASE=%.4f  OSE=%.4f\n\n', ...
    sse_cgrdpg, sse_ase, sse_ose);

%% 6. Vertex-wise coverage (all 6 methods)
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
    if mod(i, 200) == 0
        fprintf('  Vertex %d/%d\n', i, n);
    end

    % cgrdpg
    err_cgrdpg = X0(i,:) - X_cgrdpg(i,:);
    G_true_cgrdpg = compute_fisher_info_cgrdpg(i, X0, Y0, Z0, tau);
    G_plugin_cgrdpg = compute_fisher_info_cgrdpg(i, X_cgrdpg, Y_cgrdpg, Z_cgrdpg, tau);
    results.cgrdpg_true(i) = check_coverage(err_cgrdpg, G_true_cgrdpg, chi2_crit, n + p_cov);
    results.cgrdpg_plugin(i) = check_coverage(err_cgrdpg, G_plugin_cgrdpg, chi2_crit, n + p_cov);

    % ASE
    err_ase = X_ase(i,:) - X0(i,:);
    Prec_true_ase = compute_prec_ase(i, X0, S, eps_clip);
    Prec_plugin_ase = compute_prec_ase(i, X_ase, S_estimated, eps_clip);
    results.ase_true(i) = check_coverage(err_ase, Prec_true_ase, chi2_crit, 1.0);
    results.ase_plugin(i) = check_coverage(err_ase, Prec_plugin_ase, chi2_crit, 1.0);

    % OSE
    err_ose = X_ose(i,:) - X0(i,:);
    Prec_true_ose = compute_prec_ose(i, X0, Y0, eps_clip);
    Prec_plugin_ose = compute_prec_ose(i, X_ose, X_ose * S_estimated, eps_clip);
    results.ose_true(i) = check_coverage(err_ose, Prec_true_ose, chi2_crit, 1.0);
    results.ose_plugin(i) = check_coverage(err_ose, Prec_plugin_ose, chi2_crit, 1.0);
end

cov_time = toc(t0);

% Overall coverage rates
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

%% 7. Save results
out_file = fullfile(output_dir, sprintf('rep_%03d.mat', rep_id));
save(out_file, ...
    'rep_id', 'n', 'p_cov', 'd', 'tau', 'S', 'S_estimated', ...
    'results', 'overall_cov', 'n_na', ...
    'sse_cgrdpg', 'sse_ase', 'sse_ose', ...
    'cgrdpg_time', 'ase_time', 'ose_time', 'cov_time', 'rep_time', ...
    'exitflag', 'output', 'step_size_init');

fprintf('\nResults saved to: %s\n', out_file);
fprintf('============================================================================\n');

end
