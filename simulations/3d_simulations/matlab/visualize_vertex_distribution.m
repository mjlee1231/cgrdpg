% Visualize empirical vs theoretical distribution for a single vertex
% Compare estimated positions across 100 reps with theoretical normal density
clear; clc;
addpath('core');

fprintf('Visualizing cgrdpg Vertex Distribution: Empirical vs Theoretical\n');
fprintf('================================================================\n\n');

% Parameters
results_dir = 'results_matlab_3d_coverage';  % Batch surrogate results
n_reps = 100;
vertex_id = randi(1000);  % Random vertex, or set specific: vertex_id = 500;

fprintf('Selected vertex: %d\n', vertex_id);
fprintf('Loading %d replications...\n\n', n_reps);

% Call the visualization function
visualize_single_vertex(vertex_id, results_dir, n_reps);

fprintf('\nVisualization complete for vertex %d\n', vertex_id);
fprintf('Plots saved in: vertex_distribution_plots/\n');

fprintf('\nDiagnosis:\n');
fprintf('  If empirical distribution is WIDER than theoretical ellipse:\n');
fprintf('    → Fisher information underestimates variance (too optimistic)\n');
fprintf('    → Results in lower coverage than 95%%\n\n');
fprintf('  If empirical distribution is NARROWER than theoretical ellipse:\n');
fprintf('    → Fisher information overestimates variance (too conservative)\n');
fprintf('    → Results in higher coverage than 95%%\n\n');
