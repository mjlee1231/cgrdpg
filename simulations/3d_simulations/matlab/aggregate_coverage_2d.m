% Aggregate coverage results from 100 replications - VERTEX-WISE
clear; clc;

fprintf('Aggregating 2D Coverage Results (Vertex-wise)\n');
fprintf('========================================\n\n');

results_dir = 'results_matlab_2d_coverage';
n_reps = 100;
n_vertices = 500;  % From simulation parameters

% Initialize storage: (n_vertices x n_reps)
coverage_true_mat = nan(n_vertices, n_reps);
coverage_plugin_mat = nan(n_vertices, n_reps);
all_sse = nan(n_reps, 1);
all_times = nan(n_reps, 1);
all_converged = nan(n_reps, 1);

% Load each replication
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);

        % Store vertex-wise coverage indicators (n_vertices x 1)
        coverage_true_mat(:, rep) = data.results.cgrdpg_true;
        coverage_plugin_mat(:, rep) = data.results.cgrdpg_plugin;

        all_sse(rep) = data.sse_cgrdpg;
        all_times(rep) = data.cgrdpg_time;
        all_converged(rep) = (data.exitflag == 1);
        n_loaded = n_loaded + 1;
    else
        fprintf('Warning: Rep %d not found\n', rep);
    end
end

fprintf('Loaded %d/%d replications\n\n', n_loaded, n_reps);

%% Compute VERTEX-WISE coverage rates (across replications)
% For each vertex: mean across 100 reps
vertexwise_coverage_true = mean(coverage_true_mat, 2, 'omitnan');    % (n_vertices x 1)
vertexwise_coverage_plugin = mean(coverage_plugin_mat, 2, 'omitnan'); % (n_vertices x 1)

% Overall coverage (average across all vertices)
overall_coverage_true = mean(vertexwise_coverage_true, 'omitnan');
overall_coverage_plugin = mean(vertexwise_coverage_plugin, 'omitnan');

%% Summary statistics
fprintf('========================================\n');
fprintf('Vertex-wise Coverage Summary\n');
fprintf('========================================\n\n');

fprintf('cgrdpg-TRUE:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_true);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_true), 100 * max(vertexwise_coverage_true));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_true, 'omitnan'));

fprintf('cgrdpg-PLUGIN:\n');
fprintf('  Overall (avg across vertices): %.2f%%\n', 100 * overall_coverage_plugin);
fprintf('  Vertex coverage range:         [%.2f%%, %.2f%%]\n', ...
    100 * min(vertexwise_coverage_plugin), 100 * max(vertexwise_coverage_plugin));
fprintf('  Std across vertices:           %.2f%%\n\n', 100 * std(vertexwise_coverage_plugin, 'omitnan'));

fprintf('Optimization:\n');
fprintf('  Mean SSE:        %.4f\n', mean(all_sse, 'omitnan'));
fprintf('  Convergence:     %.1f%%\n', 100 * mean(all_converged, 'omitnan'));
fprintf('  Mean time:       %.1f sec\n', mean(all_times, 'omitnan'));

%% Save aggregated results
save(fullfile(results_dir, 'aggregated_vertexwise_results.mat'), ...
    'vertexwise_coverage_true', 'vertexwise_coverage_plugin', ...
    'coverage_true_mat', 'coverage_plugin_mat', ...
    'overall_coverage_true', 'overall_coverage_plugin', ...
    'all_sse', 'all_times', 'all_converged', ...
    'n_loaded', 'n_reps', 'n_vertices');

fprintf('\nVertex-wise results saved to: %s\n', ...
    fullfile(results_dir, 'aggregated_vertexwise_results.mat'));
