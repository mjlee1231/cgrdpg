% Aggregate coverage results from 100 replications
clear; clc;

fprintf('Aggregating 2D Coverage Results\n');
fprintf('========================================\n\n');

results_dir = 'results_matlab_2d_coverage';
n_reps = 100;

% Initialize storage
all_coverage_true = nan(n_reps, 1);
all_coverage_plugin = nan(n_reps, 1);
all_sse = nan(n_reps, 1);
all_times = nan(n_reps, 1);
all_converged = nan(n_reps, 1);

% Load each replication
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);
        all_coverage_true(rep) = data.overall_cov.cgrdpg_true;
        all_coverage_plugin(rep) = data.overall_cov.cgrdpg_plugin;
        all_sse(rep) = data.sse_cgrdpg;
        all_times(rep) = data.cgrdpg_time;
        all_converged(rep) = (data.exitflag == 1);
        n_loaded = n_loaded + 1;
    else
        fprintf('Warning: Rep %d not found\n', rep);
    end
end

fprintf('Loaded %d/%d replications\n\n', n_loaded, n_reps);

%% Summary statistics
fprintf('========================================\n');
fprintf('Coverage Summary (across %d reps)\n', n_loaded);
fprintf('========================================\n\n');

fprintf('cgrdpg-TRUE:\n');
fprintf('  Mean coverage: %.2f%%\n', 100 * mean(all_coverage_true, 'omitnan'));
fprintf('  Std:           %.2f%%\n', 100 * std(all_coverage_true, 'omitnan'));
fprintf('  Min:           %.2f%%\n', 100 * min(all_coverage_true));
fprintf('  Max:           %.2f%%\n\n', 100 * max(all_coverage_true));

fprintf('cgrdpg-PLUGIN:\n');
fprintf('  Mean coverage: %.2f%%\n', 100 * mean(all_coverage_plugin, 'omitnan'));
fprintf('  Std:           %.2f%%\n', 100 * std(all_coverage_plugin, 'omitnan'));
fprintf('  Min:           %.2f%%\n', 100 * min(all_coverage_plugin));
fprintf('  Max:           %.2f%%\n\n', 100 * max(all_coverage_plugin));

fprintf('Optimization:\n');
fprintf('  Mean SSE:        %.4f\n', mean(all_sse, 'omitnan'));
fprintf('  Convergence:     %.1f%%\n', 100 * mean(all_converged, 'omitnan'));
fprintf('  Mean time:       %.1f sec\n', mean(all_times, 'omitnan'));

%% Save aggregated results
save(fullfile(results_dir, 'aggregated_results.mat'), ...
    'all_coverage_true', 'all_coverage_plugin', 'all_sse', ...
    'all_times', 'all_converged', 'n_loaded', 'n_reps');

fprintf('\nAggregated results saved to: %s\n', ...
    fullfile(results_dir, 'aggregated_results.mat'));
