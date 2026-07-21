% Visualize specific vertices of interest
% Edit the vertex_ids list below to specify which vertices to plot
clear; clc;
addpath('core');

fprintf('Visualizing cgrdpg Vertex Distribution: Specific Vertices\n');
fprintf('=========================================================\n\n');

% Parameters
results_dir = 'results_matlab_3d_coverage';
n_reps = 100;

% EDIT THIS LINE to specify which vertices to visualize
vertex_ids = [1, 100, 250, 500, 750, 1000];  % Example: first, middle, last, etc.

fprintf('Creating plots for %d vertices: [%s]\n\n', ...
    length(vertex_ids), strjoin(string(vertex_ids), ', '));

% Create plots for each vertex
for v_idx = 1:length(vertex_ids)
    vertex_id = vertex_ids(v_idx);
    fprintf('Processing vertex %d (%d/%d)...\n', vertex_id, v_idx, length(vertex_ids));

    try
        visualize_single_vertex(vertex_id, results_dir, n_reps);
        fprintf('  ✓ Completed vertex %d\n\n', vertex_id);
    catch ME
        fprintf('  ✗ Error for vertex %d: %s\n\n', vertex_id, ME.message);
    end
end

fprintf('All visualizations complete!\n');
fprintf('Plots saved in: vertex_distribution_plots/\n');
