% Aggregate ASE/OSE/cgrdpg coverage results from 100 replications (n=2000)
clear; clc;

fprintf('Aggregating 3D ASE/OSE/cgrdpg Coverage Results (n=2000, p_cov=1000)\n');
fprintf('====================================================================\n\n');

results_dir = 'results_matlab_3d_coverage_n2000';
n_reps = 100;
n_vertices = 2000;
d = 3;

% Compute population eigendecomposition for alignment analysis
t = (1:n_vertices)' / n_vertices;
X0 = [0.42*t + 0.46, ...
      0.27 * sin(2*pi*t) + 0.46, ...
      0.20 * cos(4*pi*t)];
S = diag([1, 1, -1]);
P_pop = X0 * S * X0';

% Population eigendecomposition (sorted by magnitude)
[V_pop, D_pop] = eig(P_pop);
eigvals_pop = diag(D_pop);
[~, idx_pop] = sort(abs(eigvals_pop), 'descend');
eigvals_pop_sorted = eigvals_pop(idx_pop);
V_pop_sorted = V_pop(:, idx_pop);

% Population 3rd eigenvector (for alignment analysis)
v3_pop = V_pop_sorted(:, 3);
fprintf('Population eigenvalues (sorted by magnitude):\n');
fprintf('  λ1 = %+.2f, λ2 = %+.2f, λ3 = %+.2f\n\n', ...
    eigvals_pop_sorted(1), eigvals_pop_sorted(2), eigvals_pop_sorted(3));

% Initialize storage: (n_vertices x n_reps)
coverage_cgrdpg_true_mat = nan(n_vertices, n_reps);
coverage_cgrdpg_plugin_mat = nan(n_vertices, n_reps);
coverage_ase_true_mat = nan(n_vertices, n_reps);
coverage_ase_plugin_mat = nan(n_vertices, n_reps);
coverage_ose_true_mat = nan(n_vertices, n_reps);
coverage_ose_plugin_mat = nan(n_vertices, n_reps);

all_sse_cgrdpg = nan(n_reps, 1);
all_sse_ase = nan(n_reps, 1);
all_sse_ose = nan(n_reps, 1);
all_times_cgrdpg = nan(n_reps, 1);
all_times_ase = nan(n_reps, 1);
all_times_ose = nan(n_reps, 1);
all_converged = nan(n_reps, 1);

% Storage for ASE eigenvalue analysis
all_ase_eigenvalues = nan(n_reps, 3);  % Top 3 eigenvalues
all_ase_eigenvector_alignments = nan(n_reps, 1);  % Alignment with population v3

% Load each replication
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);

        % Store vertex-wise coverage indicators (n_vertices x 1)
        coverage_cgrdpg_true_mat(:, rep) = data.results.cgrdpg_true;
        coverage_cgrdpg_plugin_mat(:, rep) = data.results.cgrdpg_plugin;
        coverage_ase_true_mat(:, rep) = data.results.ase_true;
        coverage_ase_plugin_mat(:, rep) = data.results.ase_plugin;
        coverage_ose_true_mat(:, rep) = data.results.ose_true;
        coverage_ose_plugin_mat(:, rep) = data.results.ose_plugin;

        all_sse_cgrdpg(rep) = data.sse_cgrdpg;
        all_sse_ase(rep) = data.sse_ase;
        all_sse_ose(rep) = data.sse_ose;
        all_times_cgrdpg(rep) = data.cgrdpg_time;
        all_times_ase(rep) = data.ase_time;
        all_times_ose(rep) = data.ose_time;
        all_converged(rep) = (data.exitflag == 1);

        % ASE eigenvalue analysis
        if isfield(data, 'ase_eigenvalues') && isfield(data, 'ase_eigenvectors')
            all_ase_eigenvalues(rep, :) = data.ase_eigenvalues';
            % Compute alignment with population 3rd eigenvector
            v3_obs = data.ase_eigenvectors(:, 3);
            all_ase_eigenvector_alignments(rep) = abs(v3_pop' * v3_obs);
        end

        n_loaded = n_loaded + 1;
    else
        fprintf('Warning: Rep %d not found\n', rep);
    end
end

fprintf('Loaded %d/%d replications\n\n', n_loaded, n_reps);

%% Compute VERTEX-WISE coverage rates (across replications)
% For each vertex: mean across 100 reps
vertexwise_coverage_cgrdpg_true = mean(coverage_cgrdpg_true_mat, 2, 'omitnan');
vertexwise_coverage_cgrdpg_plugin = mean(coverage_cgrdpg_plugin_mat, 2, 'omitnan');
vertexwise_coverage_ase_true = mean(coverage_ase_true_mat, 2, 'omitnan');
vertexwise_coverage_ase_plugin = mean(coverage_ase_plugin_mat, 2, 'omitnan');
vertexwise_coverage_ose_true = mean(coverage_ose_true_mat, 2, 'omitnan');
vertexwise_coverage_ose_plugin = mean(coverage_ose_plugin_mat, 2, 'omitnan');

% Overall coverage (average across all vertices)
overall_coverage_cgrdpg_true = mean(vertexwise_coverage_cgrdpg_true, 'omitnan');
overall_coverage_cgrdpg_plugin = mean(vertexwise_coverage_cgrdpg_plugin, 'omitnan');
overall_coverage_ase_true = mean(vertexwise_coverage_ase_true, 'omitnan');
overall_coverage_ase_plugin = mean(vertexwise_coverage_ase_plugin, 'omitnan');
overall_coverage_ose_true = mean(vertexwise_coverage_ose_true, 'omitnan');
overall_coverage_ose_plugin = mean(vertexwise_coverage_ose_plugin, 'omitnan');

%% Summary statistics
fprintf('========================================\n');
fprintf('Vertex-wise Coverage Summary (3D, n=2000)\n');
fprintf('========================================\n\n');

fprintf('cgrdpg-TRUE:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_cgrdpg_true);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_cgrdpg_true), 100 * max(vertexwise_coverage_cgrdpg_true));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_cgrdpg_true, 'omitnan'));

fprintf('cgrdpg-PLUGIN:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_cgrdpg_plugin);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_cgrdpg_plugin), 100 * max(vertexwise_coverage_cgrdpg_plugin));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_cgrdpg_plugin, 'omitnan'));

fprintf('ASE-TRUE:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_ase_true);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_ase_true), 100 * max(vertexwise_coverage_ase_true));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_ase_true, 'omitnan'));

fprintf('ASE-PLUGIN:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_ase_plugin);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_ase_plugin), 100 * max(vertexwise_coverage_ase_plugin));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_ase_plugin, 'omitnan'));

fprintf('OSE-TRUE:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_ose_true);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_ose_true), 100 * max(vertexwise_coverage_ose_true));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_ose_true, 'omitnan'));

fprintf('OSE-PLUGIN:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_ose_plugin);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_ose_plugin), 100 * max(vertexwise_coverage_ose_plugin));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_ose_plugin, 'omitnan'));

fprintf('Optimization:\n');
fprintf('  cgrdpg Convergence: %.1f%%\n', 100 * mean(all_converged, 'omitnan'));
fprintf('  cgrdpg Mean time:   %.1f (SD: %.1f) sec\n', mean(all_times_cgrdpg, 'omitnan'), std(all_times_cgrdpg, 'omitnan'));
fprintf('  ASE Mean time:      %.1f (SD: %.1f) sec\n', mean(all_times_ase, 'omitnan'), std(all_times_ase, 'omitnan'));
fprintf('  OSE Mean time:      %.1f (SD: %.1f) sec\n\n', mean(all_times_ose, 'omitnan'), std(all_times_ose, 'omitnan'));

fprintf('SSE Distribution (cgrdpg):\n');
fprintf('  Mean:            %.4f\n', mean(all_sse_cgrdpg, 'omitnan'));
fprintf('  Std:             %.4f\n', std(all_sse_cgrdpg, 'omitnan'));
fprintf('  Min:             %.4f\n', min(all_sse_cgrdpg));
fprintf('  25th percentile: %.4f\n', quantile(all_sse_cgrdpg, 0.25));
fprintf('  Median:          %.4f\n', median(all_sse_cgrdpg, 'omitnan'));
fprintf('  75th percentile: %.4f\n', quantile(all_sse_cgrdpg, 0.75));
fprintf('  Max:             %.4f\n\n', max(all_sse_cgrdpg));

fprintf('SSE Distribution (ASE):\n');
fprintf('  Mean:            %.4f\n', mean(all_sse_ase, 'omitnan'));
fprintf('  Std:             %.4f\n', std(all_sse_ase, 'omitnan'));
fprintf('  Min:             %.4f\n', min(all_sse_ase));
fprintf('  25th percentile: %.4f\n', quantile(all_sse_ase, 0.25));
fprintf('  Median:          %.4f\n', median(all_sse_ase, 'omitnan'));
fprintf('  75th percentile: %.4f\n', quantile(all_sse_ase, 0.75));
fprintf('  Max:             %.4f\n\n', max(all_sse_ase));

fprintf('SSE Distribution (OSE):\n');
fprintf('  Mean:            %.4f\n', mean(all_sse_ose, 'omitnan'));
fprintf('  Std:             %.4f\n', std(all_sse_ose, 'omitnan'));
fprintf('  Min:             %.4f\n', min(all_sse_ose));
fprintf('  25th percentile: %.4f\n', quantile(all_sse_ose, 0.25));
fprintf('  Median:          %.4f\n', median(all_sse_ose, 'omitnan'));
fprintf('  75th percentile: %.4f\n', quantile(all_sse_ose, 0.75));
fprintf('  Max:             %.4f\n', max(all_sse_ose));

%% ASE Eigenvalue Analysis
fprintf('\n========================================\n');
fprintf('ASE Eigenvalue Analysis (Indefinite GRDPG)\n');
fprintf('========================================\n\n');

% Count sign patterns of observed rank 3 eigenvalue
n_eigenvalue_available = sum(~isnan(all_ase_eigenvalues(:, 3)));
signs_rank3 = sign(all_ase_eigenvalues(:, 3));
n_positive = sum(signs_rank3 > 0);
n_negative = sum(signs_rank3 < 0);
n_zero = sum(signs_rank3 == 0);

fprintf('Observed Rank 3 Eigenvalue Signs (%d reps with data):\n', n_eigenvalue_available);
fprintf('  Positive: %3d (%.1f%%)\n', n_positive, 100 * n_positive / n_eigenvalue_available);
fprintf('  Negative: %3d (%.1f%%)\n', n_negative, 100 * n_negative / n_eigenvalue_available);
fprintf('  Zero:     %3d (%.1f%%)\n\n', n_zero, 100 * n_zero / n_eigenvalue_available);

fprintf('Expected: Negative (matches population λ3 = -20.00)\n\n');

% Eigenvalue magnitude statistics
fprintf('Observed Eigenvalue Magnitudes (across %d reps):\n', n_eigenvalue_available);
for k = 1:3
    vals = abs(all_ase_eigenvalues(:, k));
    fprintf('  |λ%d|: Mean=%.2f, SD=%.2f, Range=[%.2f, %.2f]\n', ...
        k, mean(vals, 'omitnan'), std(vals, 'omitnan'), min(vals), max(vals));
end
fprintf('\n');

% Eigenvector alignment analysis
n_alignment_available = sum(~isnan(all_ase_eigenvector_alignments));
fprintf('Eigenvector Alignment with Population v3 (%d reps with data):\n', n_alignment_available);
fprintf('  Mean alignment: %.4f\n', mean(all_ase_eigenvector_alignments, 'omitnan'));
fprintf('  SD:             %.4f\n', std(all_ase_eigenvector_alignments, 'omitnan'));
fprintf('  Min:            %.4f\n', min(all_ase_eigenvector_alignments));
fprintf('  Max:            %.4f\n\n', max(all_ase_eigenvector_alignments));

% Identify problematic cases
idx_positive_rank3 = find(signs_rank3 > 0);
if ~isempty(idx_positive_rank3)
    fprintf('⚠️  WARNING: %d replications have POSITIVE rank 3 eigenvalue\n', length(idx_positive_rank3));
    fprintf('   Reps with positive λ3: %s\n', mat2str(idx_positive_rank3'));
    fprintf('   Mean alignment for these reps: %.4f\n', ...
        mean(all_ase_eigenvector_alignments(idx_positive_rank3), 'omitnan'));
    fprintf('   (Low alignment suggests wrong eigenvector selection)\n\n');
else
    fprintf('✓ All replications have NEGATIVE rank 3 eigenvalue (correct sign)\n\n');
end

%% Save aggregated results
save(fullfile(results_dir, 'aggregated_vertexwise_results.mat'), ...
    'vertexwise_coverage_cgrdpg_true', 'vertexwise_coverage_cgrdpg_plugin', ...
    'vertexwise_coverage_ase_true', 'vertexwise_coverage_ase_plugin', ...
    'vertexwise_coverage_ose_true', 'vertexwise_coverage_ose_plugin', ...
    'coverage_cgrdpg_true_mat', 'coverage_cgrdpg_plugin_mat', ...
    'coverage_ase_true_mat', 'coverage_ase_plugin_mat', ...
    'coverage_ose_true_mat', 'coverage_ose_plugin_mat', ...
    'overall_coverage_cgrdpg_true', 'overall_coverage_cgrdpg_plugin', ...
    'overall_coverage_ase_true', 'overall_coverage_ase_plugin', ...
    'overall_coverage_ose_true', 'overall_coverage_ose_plugin', ...
    'all_sse_cgrdpg', 'all_sse_ase', 'all_sse_ose', ...
    'all_times_cgrdpg', 'all_times_ase', 'all_times_ose', ...
    'all_converged', 'n_loaded', 'n_reps', 'n_vertices', ...
    'all_ase_eigenvalues', 'all_ase_eigenvector_alignments', ...
    'eigvals_pop_sorted', 'v3_pop');

fprintf('\nVertex-wise results saved to: %s\n', ...
    fullfile(results_dir, 'aggregated_vertexwise_results.mat'));
