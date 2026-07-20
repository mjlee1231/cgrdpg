function coverage_3d_single_rep_one_step(rep_id)
% COVERAGE_3D_SINGLE_REP_ONE_STEP Test one-step cgrdpg estimators
%
% Compares:
%   - ASE (baseline)
%   - cgrdpg-ONE-STEP-BATCH (batch Newton-Raphson, 1 iteration)
%   - cgrdpg-ONE-STEP-JACOBI (Jacobi Newton-Raphson, 1 iteration)
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
tau = 0.001;
eps_clip = 1e-10;
chi2_crit = chi2inv(0.95, d);
step_size_jacobi = 1.0;  % Full Newton step for one-step

% Output directory
output_dir = 'results_matlab_3d_one_step';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('============================================================================\n');
fprintf('  One-Step cgrdpg Estimators: 3D GRDPG, n=%d\n', n);
fprintf('  Replication %d/100\n', rep_id);
fprintf('  Methods: ASE, cgrdpg-ONE-STEP-BATCH, cgrdpg-ONE-STEP-JACOBI\n');
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

%% 3. ASE (initialization for one-step methods)
fprintf('Computing ASE...\n');
t0 = tic;
[X_ase_unsigned, X_ase_signed, S_estimated] = fit_ase(A, d, p);
ase_time = toc(t0);
[X_ase, ~] = procrustes_align(X_ase_unsigned, X0);
fprintf('ASE: time=%.1fs, S_est=diag([%+d,%+d,%+d])\n', ...
    ase_time, diag(S_estimated));

%% 4. cgrdpg-ONE-STEP-BATCH
fprintf('Computing cgrdpg-ONE-STEP-BATCH...\n');
t0 = tic;
[X_one_step_batch, Z_one_step_batch] = ...
    compute_cgrdpg_one_step_batch(A, B, X_ase_unsigned, S_estimated, tau);
one_step_batch_time = toc(t0);
[X_one_step_batch_aligned, ~] = procrustes_align(X_one_step_batch, X0);
Y_one_step_batch = X_one_step_batch_aligned * S_estimated;
fprintf('cgrdpg-ONE-STEP-BATCH: time=%.1fs\n', one_step_batch_time);

%% 5. cgrdpg-ONE-STEP-JACOBI
fprintf('Computing cgrdpg-ONE-STEP-JACOBI...\n');
t0 = tic;
[X_one_step_jacobi, Z_one_step_jacobi] = ...
    compute_cgrdpg_one_step_jacobi(A, B, X_ase_unsigned, S_estimated, tau, step_size_jacobi);
one_step_jacobi_time = toc(t0);
[X_one_step_jacobi_aligned, ~] = procrustes_align(X_one_step_jacobi, X0);
Y_one_step_jacobi = X_one_step_jacobi_aligned * S_estimated;
fprintf('cgrdpg-ONE-STEP-JACOBI: time=%.1fs\n', one_step_jacobi_time);

%% 6. SSE
sse_ase = sum((X_ase - X0).^2, 'all');
sse_one_step_batch = sum((X_one_step_batch_aligned - X0).^2, 'all');
sse_one_step_jacobi = sum((X_one_step_jacobi_aligned - X0).^2, 'all');

fprintf('\nSSE:\n');
fprintf('  ASE:                %.4f (baseline)\n', sse_ase);
fprintf('  cgrdpg-ONE-STEP-BATCH:  %.4f', sse_one_step_batch);
if sse_one_step_batch < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n', 100*(sse_ase - sse_one_step_batch)/sse_ase);
else
    fprintf(' (%.1f%% worse)\n', 100*(sse_one_step_batch - sse_ase)/sse_ase);
end
fprintf('  cgrdpg-ONE-STEP-JACOBI: %.4f', sse_one_step_jacobi);
if sse_one_step_jacobi < sse_ase
    fprintf(' ✓ (%.1f%% improvement)\n\n', 100*(sse_ase - sse_one_step_jacobi)/sse_ase);
else
    fprintf(' (%.1f%% worse)\n\n', 100*(sse_one_step_jacobi - sse_ase)/sse_ase);
end

%% 7. Vertex-wise coverage (all methods with TRUE and PLUGIN)
fprintf('Computing vertex-wise coverage for all %d vertices...\n', n);

results = struct();
results.ase_true = nan(n, 1);
results.ase_plugin = nan(n, 1);
results.one_step_batch_true = nan(n, 1);
results.one_step_batch_plugin = nan(n, 1);
results.one_step_jacobi_true = nan(n, 1);
results.one_step_jacobi_plugin = nan(n, 1);

t0 = tic;
for i = 1:n
    if mod(i, 200) == 0
        fprintf('  Vertex %d/%d\n', i, n);
    end

    % ASE
    err_ase = X_ase(i,:) - X0(i,:);
    Prec_true_ase = compute_prec_ase(i, X0, S, eps_clip);
    Prec_plugin_ase = compute_prec_ase(i, X_ase, S_estimated, eps_clip);
    results.ase_true(i) = check_coverage(err_ase, Prec_true_ase, chi2_crit, 1.0);
    results.ase_plugin(i) = check_coverage(err_ase, Prec_plugin_ase, chi2_crit, 1.0);

    % cgrdpg-ONE-STEP-BATCH
    err_batch = X0(i,:) - X_one_step_batch_aligned(i,:);
    G_true_batch = compute_fisher_info_cgrdpg(i, X0, Y0, Z0, tau);
    G_plugin_batch = compute_fisher_info_cgrdpg(i, X_one_step_batch_aligned, ...
        Y_one_step_batch, Z_one_step_batch, tau);
    results.one_step_batch_true(i) = check_coverage(err_batch, G_true_batch, chi2_crit, n + p_cov);
    results.one_step_batch_plugin(i) = check_coverage(err_batch, G_plugin_batch, chi2_crit, n + p_cov);

    % cgrdpg-ONE-STEP-JACOBI
    err_jacobi = X0(i,:) - X_one_step_jacobi_aligned(i,:);
    G_true_jacobi = compute_fisher_info_cgrdpg(i, X0, Y0, Z0, tau);
    G_plugin_jacobi = compute_fisher_info_cgrdpg(i, X_one_step_jacobi_aligned, ...
        Y_one_step_jacobi, Z_one_step_jacobi, tau);
    results.one_step_jacobi_true(i) = check_coverage(err_jacobi, G_true_jacobi, chi2_crit, n + p_cov);
    results.one_step_jacobi_plugin(i) = check_coverage(err_jacobi, G_plugin_jacobi, chi2_crit, n + p_cov);
end

cov_time = toc(t0);

% Overall coverage rates
overall_cov = struct();
overall_cov.ase_true = mean(results.ase_true, 'omitnan');
overall_cov.ase_plugin = mean(results.ase_plugin, 'omitnan');
overall_cov.one_step_batch_true = mean(results.one_step_batch_true, 'omitnan');
overall_cov.one_step_batch_plugin = mean(results.one_step_batch_plugin, 'omitnan');
overall_cov.one_step_jacobi_true = mean(results.one_step_jacobi_true, 'omitnan');
overall_cov.one_step_jacobi_plugin = mean(results.one_step_jacobi_plugin, 'omitnan');

% Count NAs
n_na = struct();
n_na.ase_true = sum(isnan(results.ase_true));
n_na.ase_plugin = sum(isnan(results.ase_plugin));
n_na.one_step_batch_true = sum(isnan(results.one_step_batch_true));
n_na.one_step_batch_plugin = sum(isnan(results.one_step_batch_plugin));
n_na.one_step_jacobi_true = sum(isnan(results.one_step_jacobi_true));
n_na.one_step_jacobi_plugin = sum(isnan(results.one_step_jacobi_plugin));

rep_time = toc(rep_start) / 60;  % minutes

fprintf('\nOverall coverage (this rep):\n');
fprintf('  ASE-TRUE                    %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.ase_true, n_na.ase_true);
fprintf('  ASE-PLUGIN                  %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.ase_plugin, n_na.ase_plugin);
fprintf('  cgrdpg-ONE-STEP-BATCH-TRUE      %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.one_step_batch_true, n_na.one_step_batch_true);
fprintf('  cgrdpg-ONE-STEP-BATCH-PLUGIN    %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.one_step_batch_plugin, n_na.one_step_batch_plugin);
fprintf('  cgrdpg-ONE-STEP-JACOBI-TRUE     %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.one_step_jacobi_true, n_na.one_step_jacobi_true);
fprintf('  cgrdpg-ONE-STEP-JACOBI-PLUGIN   %.1f%%  (NAs: %d)\n', ...
    100 * overall_cov.one_step_jacobi_plugin, n_na.one_step_jacobi_plugin);
fprintf('\nTotal rep time: %.2f min\n', rep_time);

%% 8. Save results
out_file = fullfile(output_dir, sprintf('rep_%03d.mat', rep_id));
save(out_file, ...
    'rep_id', 'n', 'p_cov', 'd', 'tau', 'S', 'S_estimated', ...
    'results', 'overall_cov', 'n_na', ...
    'sse_ase', 'sse_one_step_batch', 'sse_one_step_jacobi', ...
    'ase_time', 'one_step_batch_time', 'one_step_jacobi_time', ...
    'cov_time', 'rep_time', 'step_size_jacobi');

fprintf('\nResults saved to: %s\n', out_file);
fprintf('============================================================================\n');

end
