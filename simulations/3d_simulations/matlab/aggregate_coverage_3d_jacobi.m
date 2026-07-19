% Aggregate Jacobi method coverage results from 100 replications
clear; clc;

fprintf('Aggregating 3D Coverage Results (JACOBI METHOD)\n');
fprintf('===============================================\n\n');

results_dir = 'results_matlab_3d_coverage_jacobi';
n_reps = 100;
n_vertices = 1000;

% Initialize storage
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

% Load each replication
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);

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
        n_loaded = n_loaded + 1;
    else
        fprintf('Warning: Rep %d not found\n', rep);
    end
end

fprintf('Loaded %d/%d replications\n\n', n_loaded, n_reps);

%% Compute vertex-wise coverage rates
vertexwise_coverage_cgrdpg_true = mean(coverage_cgrdpg_true_mat, 2, 'omitnan');
vertexwise_coverage_cgrdpg_plugin = mean(coverage_cgrdpg_plugin_mat, 2, 'omitnan');
vertexwise_coverage_ase_true = mean(coverage_ase_true_mat, 2, 'omitnan');
vertexwise_coverage_ase_plugin = mean(coverage_ase_plugin_mat, 2, 'omitnan');
vertexwise_coverage_ose_true = mean(coverage_ose_true_mat, 2, 'omitnan');
vertexwise_coverage_ose_plugin = mean(coverage_ose_plugin_mat, 2, 'omitnan');

% Overall coverage
overall_coverage_cgrdpg_true = mean(vertexwise_coverage_cgrdpg_true, 'omitnan');
overall_coverage_cgrdpg_plugin = mean(vertexwise_coverage_cgrdpg_plugin, 'omitnan');
overall_coverage_ase_true = mean(vertexwise_coverage_ase_true, 'omitnan');
overall_coverage_ase_plugin = mean(vertexwise_coverage_ase_plugin, 'omitnan');
overall_coverage_ose_true = mean(vertexwise_coverage_ose_true, 'omitnan');
overall_coverage_ose_plugin = mean(vertexwise_coverage_ose_plugin, 'omitnan');

%% Summary
fprintf('========================================\n');
fprintf('Vertex-wise Coverage Summary (JACOBI METHOD)\n');
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
fprintf('  cgrdpg Mean time:   %.1f sec\n', mean(all_times_cgrdpg, 'omitnan'));
fprintf('  ASE Mean time:      %.1f sec\n', mean(all_times_ase, 'omitnan'));
fprintf('  OSE Mean time:      %.1f sec\n\n', mean(all_times_ose, 'omitnan'));

fprintf('SSE Distribution (cgrdpg with Jacobi):\n');
fprintf('  Mean:   %.4f\n', mean(all_sse_cgrdpg, 'omitnan'));
fprintf('  Median: %.4f\n', median(all_sse_cgrdpg, 'omitnan'));
fprintf('  Std:    %.4f\n', std(all_sse_cgrdpg, 'omitnan'));
fprintf('  Min:    %.4f\n', min(all_sse_cgrdpg));
fprintf('  Max:    %.4f\n\n', max(all_sse_cgrdpg));

fprintf('========================================\n');
fprintf('Comparison with Batch Surrogate\n');
fprintf('========================================\n');
fprintf('Batch cgrdpg SSE mean:   48.60\n');
fprintf('Jacobi cgrdpg SSE mean:  %.2f\n', mean(all_sse_cgrdpg, 'omitnan'));
fprintf('Batch cgrdpg-TRUE:       26.42%%\n');
fprintf('Jacobi cgrdpg-TRUE:      %.2f%%\n', 100 * overall_coverage_cgrdpg_true);

%% Save
save(fullfile(results_dir, 'aggregated_jacobi_results.mat'), ...
    'vertexwise_coverage_cgrdpg_true', 'vertexwise_coverage_cgrdpg_plugin', ...
    'vertexwise_coverage_ase_true', 'vertexwise_coverage_ase_plugin', ...
    'vertexwise_coverage_ose_true', 'vertexwise_coverage_ose_plugin', ...
    'overall_coverage_cgrdpg_true', 'overall_coverage_cgrdpg_plugin', ...
    'overall_coverage_ase_true', 'overall_coverage_ase_plugin', ...
    'overall_coverage_ose_true', 'overall_coverage_ose_plugin', ...
    'all_sse_cgrdpg', 'all_sse_ase', 'all_sse_ose', ...
    'all_times_cgrdpg', 'all_times_ase', 'all_times_ose', ...
    'all_converged', 'n_loaded', 'n_reps', 'n_vertices');

fprintf('\nResults saved to: %s\n', ...
    fullfile(results_dir, 'aggregated_jacobi_results.mat'));
