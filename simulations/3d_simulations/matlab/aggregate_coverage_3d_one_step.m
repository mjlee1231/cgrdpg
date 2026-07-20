% Aggregate one-step cgrdpg coverage results from 100 replications
clear; clc;

fprintf('Aggregating 3D One-Step cgrdpg Coverage Results\n');
fprintf('================================================\n\n');

results_dir = 'results_matlab_3d_one_step';
n_reps = 100;
n_vertices = 1000;

% Initialize storage
coverage_ase_true_mat = nan(n_vertices, n_reps);
coverage_ase_plugin_mat = nan(n_vertices, n_reps);
coverage_one_step_batch_true_mat = nan(n_vertices, n_reps);
coverage_one_step_batch_plugin_mat = nan(n_vertices, n_reps);
coverage_one_step_jacobi_true_mat = nan(n_vertices, n_reps);
coverage_one_step_jacobi_plugin_mat = nan(n_vertices, n_reps);

all_sse_ase = nan(n_reps, 1);
all_sse_one_step_batch = nan(n_reps, 1);
all_sse_one_step_jacobi = nan(n_reps, 1);
all_times_ase = nan(n_reps, 1);
all_times_one_step_batch = nan(n_reps, 1);
all_times_one_step_jacobi = nan(n_reps, 1);

% Load each replication
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);

        coverage_ase_true_mat(:, rep) = data.results.ase_true;
        coverage_ase_plugin_mat(:, rep) = data.results.ase_plugin;
        coverage_one_step_batch_true_mat(:, rep) = data.results.one_step_batch_true;
        coverage_one_step_batch_plugin_mat(:, rep) = data.results.one_step_batch_plugin;
        coverage_one_step_jacobi_true_mat(:, rep) = data.results.one_step_jacobi_true;
        coverage_one_step_jacobi_plugin_mat(:, rep) = data.results.one_step_jacobi_plugin;

        all_sse_ase(rep) = data.sse_ase;
        all_sse_one_step_batch(rep) = data.sse_one_step_batch;
        all_sse_one_step_jacobi(rep) = data.sse_one_step_jacobi;
        all_times_ase(rep) = data.ase_time;
        all_times_one_step_batch(rep) = data.one_step_batch_time;
        all_times_one_step_jacobi(rep) = data.one_step_jacobi_time;
        n_loaded = n_loaded + 1;
    else
        fprintf('Warning: Rep %d not found\n', rep);
    end
end

fprintf('Loaded %d/%d replications\n\n', n_loaded, n_reps);

%% Compute vertex-wise coverage rates
vertexwise_coverage_ase_true = mean(coverage_ase_true_mat, 2, 'omitnan');
vertexwise_coverage_ase_plugin = mean(coverage_ase_plugin_mat, 2, 'omitnan');
vertexwise_coverage_one_step_batch_true = mean(coverage_one_step_batch_true_mat, 2, 'omitnan');
vertexwise_coverage_one_step_batch_plugin = mean(coverage_one_step_batch_plugin_mat, 2, 'omitnan');
vertexwise_coverage_one_step_jacobi_true = mean(coverage_one_step_jacobi_true_mat, 2, 'omitnan');
vertexwise_coverage_one_step_jacobi_plugin = mean(coverage_one_step_jacobi_plugin_mat, 2, 'omitnan');

% Overall coverage
overall_coverage_ase_true = mean(vertexwise_coverage_ase_true, 'omitnan');
overall_coverage_ase_plugin = mean(vertexwise_coverage_ase_plugin, 'omitnan');
overall_coverage_one_step_batch_true = mean(vertexwise_coverage_one_step_batch_true, 'omitnan');
overall_coverage_one_step_batch_plugin = mean(vertexwise_coverage_one_step_batch_plugin, 'omitnan');
overall_coverage_one_step_jacobi_true = mean(vertexwise_coverage_one_step_jacobi_true, 'omitnan');
overall_coverage_one_step_jacobi_plugin = mean(vertexwise_coverage_one_step_jacobi_plugin, 'omitnan');

%% Summary
fprintf('========================================\n');
fprintf('Vertex-wise Coverage Summary\n');
fprintf('========================================\n\n');

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

fprintf('cgrdpg-ONE-STEP-BATCH-TRUE:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_one_step_batch_true);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_one_step_batch_true), 100 * max(vertexwise_coverage_one_step_batch_true));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_one_step_batch_true, 'omitnan'));

fprintf('cgrdpg-ONE-STEP-BATCH-PLUGIN:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_one_step_batch_plugin);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_one_step_batch_plugin), 100 * max(vertexwise_coverage_one_step_batch_plugin));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_one_step_batch_plugin, 'omitnan'));

fprintf('cgrdpg-ONE-STEP-JACOBI-TRUE:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_one_step_jacobi_true);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_one_step_jacobi_true), 100 * max(vertexwise_coverage_one_step_jacobi_true));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_one_step_jacobi_true, 'omitnan'));

fprintf('cgrdpg-ONE-STEP-JACOBI-PLUGIN:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_one_step_jacobi_plugin);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_one_step_jacobi_plugin), 100 * max(vertexwise_coverage_one_step_jacobi_plugin));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_one_step_jacobi_plugin, 'omitnan'));

fprintf('Timing:\n');
fprintf('  ASE Mean time:               %.1f sec\n', mean(all_times_ase, 'omitnan'));
fprintf('  ONE-STEP-BATCH Mean time:    %.1f sec\n', mean(all_times_one_step_batch, 'omitnan'));
fprintf('  ONE-STEP-JACOBI Mean time:   %.1f sec\n\n', mean(all_times_one_step_jacobi, 'omitnan'));

fprintf('SSE Distribution:\n');
fprintf('  ASE:               Mean=%.4f, SD=%.4f, Median=%.4f\n', ...
    mean(all_sse_ase, 'omitnan'), std(all_sse_ase, 'omitnan'), median(all_sse_ase, 'omitnan'));
fprintf('  ONE-STEP-BATCH:    Mean=%.4f, SD=%.4f, Median=%.4f\n', ...
    mean(all_sse_one_step_batch, 'omitnan'), std(all_sse_one_step_batch, 'omitnan'), ...
    median(all_sse_one_step_batch, 'omitnan'));
fprintf('  ONE-STEP-JACOBI:   Mean=%.4f, SD=%.4f, Median=%.4f\n\n', ...
    mean(all_sse_one_step_jacobi, 'omitnan'), std(all_sse_one_step_jacobi, 'omitnan'), ...
    median(all_sse_one_step_jacobi, 'omitnan'));

fprintf('SSE Improvement over ASE:\n');
fprintf('  ONE-STEP-BATCH:  %.1f%%\n', ...
    100 * (mean(all_sse_ase, 'omitnan') - mean(all_sse_one_step_batch, 'omitnan')) / mean(all_sse_ase, 'omitnan'));
fprintf('  ONE-STEP-JACOBI: %.1f%%\n\n', ...
    100 * (mean(all_sse_ase, 'omitnan') - mean(all_sse_one_step_jacobi, 'omitnan')) / mean(all_sse_ase, 'omitnan'));

fprintf('========================================\n');
fprintf('Comparison: Batch vs Jacobi\n');
fprintf('========================================\n');
fprintf('Coverage TRUE:\n');
fprintf('  Batch:  %.2f%%\n', 100 * overall_coverage_one_step_batch_true);
fprintf('  Jacobi: %.2f%%\n', 100 * overall_coverage_one_step_jacobi_true);
fprintf('Coverage PLUGIN:\n');
fprintf('  Batch:  %.2f%%\n', 100 * overall_coverage_one_step_batch_plugin);
fprintf('  Jacobi: %.2f%%\n', 100 * overall_coverage_one_step_jacobi_plugin);
fprintf('SSE:\n');
fprintf('  Batch:  %.4f\n', mean(all_sse_one_step_batch, 'omitnan'));
fprintf('  Jacobi: %.4f\n', mean(all_sse_one_step_jacobi, 'omitnan'));

%% Save
save(fullfile(results_dir, 'aggregated_one_step_results.mat'), ...
    'vertexwise_coverage_ase_true', 'vertexwise_coverage_ase_plugin', ...
    'vertexwise_coverage_one_step_batch_true', 'vertexwise_coverage_one_step_batch_plugin', ...
    'vertexwise_coverage_one_step_jacobi_true', 'vertexwise_coverage_one_step_jacobi_plugin', ...
    'overall_coverage_ase_true', 'overall_coverage_ase_plugin', ...
    'overall_coverage_one_step_batch_true', 'overall_coverage_one_step_batch_plugin', ...
    'overall_coverage_one_step_jacobi_true', 'overall_coverage_one_step_jacobi_plugin', ...
    'all_sse_ase', 'all_sse_one_step_batch', 'all_sse_one_step_jacobi', ...
    'all_times_ase', 'all_times_one_step_batch', 'all_times_one_step_jacobi', ...
    'n_loaded', 'n_reps', 'n_vertices');

fprintf('\nResults saved to: %s\n', ...
    fullfile(results_dir, 'aggregated_one_step_results.mat'));
