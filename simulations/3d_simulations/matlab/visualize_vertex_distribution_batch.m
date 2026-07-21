% Visualize empirical vs theoretical distribution for MULTIPLE vertices
% Creates diagnostic plots for several vertices to check coverage behavior
clear; clc;
addpath('core');

fprintf('Visualizing cgrdpg Vertex Distribution: Multiple Vertices\n');
fprintf('=========================================================\n\n');

% Parameters
results_dir = 'results_matlab_3d_coverage';
n_reps = 100;

% Select vertices to visualize
% Option 1: Random sample
n_vertices_to_plot = 5;
rng(42);  % For reproducibility
vertex_ids = sort(randperm(1000, n_vertices_to_plot));

% Option 2: Specific vertices (uncomment to use)
% vertex_ids = [1, 250, 500, 750, 1000];  % Start, quartiles, end

% Option 3: Vertices with extreme coverage (uncomment after first run)
% Load aggregated results to find interesting vertices
% agg = load(fullfile(results_dir, 'aggregated_vertexwise_results.mat'));
% [~, idx_low] = mink(agg.vertexwise_coverage_cgrdpg_plugin, 3);   % Lowest coverage
% [~, idx_high] = maxk(agg.vertexwise_coverage_cgrdpg_plugin, 3);  % Highest coverage
% vertex_ids = sort([idx_low; idx_high]');

fprintf('Creating plots for %d vertices: [%s]\n\n', ...
    length(vertex_ids), strjoin(string(vertex_ids), ', '));

% Create one plot per vertex
for v_idx = 1:length(vertex_ids)
    vertex_id = vertex_ids(v_idx);
    fprintf('Processing vertex %d (%d/%d)...\n', vertex_id, v_idx, length(vertex_ids));

    % Call the single-vertex visualization function
    visualize_single_vertex(vertex_id, results_dir, n_reps);

    fprintf('  Completed vertex %d\n\n', vertex_id);
end

fprintf('All visualizations complete!\n');
fprintf('Plots saved in: vertex_distribution_plots/\n');
