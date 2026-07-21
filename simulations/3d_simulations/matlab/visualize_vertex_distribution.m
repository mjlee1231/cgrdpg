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
d = 3;
alpha = 0.05;
chi2_crit = chi2inv(1 - alpha, d);
tau = 0.001;  % Same as simulation

fprintf('Selected vertex: %d\n', vertex_id);
fprintf('Loading %d replications...\n', n_reps);

% Storage for estimates across reps
X_estimates = nan(n_reps, d);
X_true = nan(1, d);
fisher_info_plugin = nan(d, d);
coverage_plugin = nan(n_reps, 1);

% Load all replications
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);
        X_estimates(rep, :) = data.X_cgrdpg(vertex_id, :);

        if rep == 1
            X_true = data.X0(vertex_id, :);
            Y_plugin = data.Y_cgrdpg;
            Z_plugin = data.Z_cgrdpg;
            n_vertices = data.n;
            p_cov = data.p_cov;
            % Compute plugin Fisher info for this vertex
            fisher_info_plugin = compute_fisher_info_cgrdpg(vertex_id, ...
                data.X_cgrdpg, Y_plugin, Z_plugin, tau);
        end

        % Check coverage for this vertex using plugin Fisher info
        % Recompute Fisher info for each rep
        G_plugin_rep = compute_fisher_info_cgrdpg(vertex_id, ...
            data.X_cgrdpg, data.Y_cgrdpg, data.Z_cgrdpg, tau);
        err = X_true - data.X_cgrdpg(vertex_id, :);
        dist = err * G_plugin_rep * err' * (n_vertices + p_cov);
        coverage_plugin(rep) = (dist <= chi2_crit);

        n_loaded = n_loaded + 1;
    end
end

fprintf('Loaded %d/%d replications\n\n', n_loaded, n_reps);

% Compute empirical statistics
X_mean = mean(X_estimates, 1, 'omitnan');
X_cov_empirical = cov(X_estimates, 'omitrows');
empirical_coverage = mean(coverage_plugin, 'omitnan');

% Theoretical covariance from plugin Fisher info
% Cov = G^{-1} / (n + p_cov)
% Note: fisher_info_plugin was loaded from rep 1 above
Cov_theoretical = inv(fisher_info_plugin) / (n_vertices + p_cov);

fprintf('Vertex %d Summary:\n', vertex_id);
fprintf('  True position:     [%.4f, %.4f, %.4f]\n', X_true);
fprintf('  Empirical mean:    [%.4f, %.4f, %.4f]\n', X_mean);
fprintf('  Empirical bias:    [%.4f, %.4f, %.4f]\n', X_mean - X_true);
fprintf('  Empirical coverage: %.1f%%\n\n', 100 * empirical_coverage);

fprintf('Empirical Covariance:\n');
disp(X_cov_empirical);
fprintf('Theoretical Covariance (plugin):\n');
disp(Cov_theoretical);

fprintf('Covariance Ratio (Empirical / Theoretical):\n');
disp(X_cov_empirical ./ Cov_theoretical);

%% Create visualization
figure('Position', [100, 100, 1400, 900]);

% Plot 1: Dimension 1 vs 2
subplot(2, 3, 1);
hold on;
scatter(X_estimates(:, 1), X_estimates(:, 2), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
plot(X_true(1), X_true(2), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
plot(X_mean(1), X_mean(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);

% Overlay theoretical 95% confidence ellipse
theta = linspace(0, 2*pi, 100);
ellipse_scale = sqrt(chi2_crit);
Cov_12 = Cov_theoretical([1,2], [1,2]);
[V, D] = eig(Cov_12);
ellipse_pts = ellipse_scale * V * sqrt(D) * [cos(theta); sin(theta)] + X_true([1,2])';
plot(ellipse_pts(1,:), ellipse_pts(2,:), 'r-', 'LineWidth', 2);

xlabel('Dimension 1');
ylabel('Dimension 2');
title(sprintf('Vertex %d: Dim 1 vs 2', vertex_id));
legend('Estimates (100 reps)', 'True position', 'Empirical mean', ...
    '95% ellipse (plugin)', 'Location', 'best');
grid on;
hold off;

% Plot 2: Dimension 1 vs 3
subplot(2, 3, 2);
hold on;
scatter(X_estimates(:, 1), X_estimates(:, 3), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
plot(X_true(1), X_true(3), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
plot(X_mean(1), X_mean(3), 'go', 'MarkerSize', 10, 'LineWidth', 2);

% Overlay theoretical 95% confidence ellipse
Cov_13 = Cov_theoretical([1,3], [1,3]);
[V, D] = eig(Cov_13);
ellipse_pts = ellipse_scale * V * sqrt(D) * [cos(theta); sin(theta)] + X_true([1,3])';
plot(ellipse_pts(1,:), ellipse_pts(2,:), 'r-', 'LineWidth', 2);

xlabel('Dimension 1');
ylabel('Dimension 3');
title(sprintf('Vertex %d: Dim 1 vs 3', vertex_id));
legend('Estimates (100 reps)', 'True position', 'Empirical mean', ...
    '95% ellipse (plugin)', 'Location', 'best');
grid on;
hold off;

% Plot 3: Dimension 2 vs 3
subplot(2, 3, 3);
hold on;
scatter(X_estimates(:, 2), X_estimates(:, 3), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
plot(X_true(2), X_true(3), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
plot(X_mean(2), X_mean(3), 'go', 'MarkerSize', 10, 'LineWidth', 2);

% Overlay theoretical 95% confidence ellipse
Cov_23 = Cov_theoretical([2,3], [2,3]);
[V, D] = eig(Cov_23);
ellipse_pts = ellipse_scale * V * sqrt(D) * [cos(theta); sin(theta)] + X_true([2,3])';
plot(ellipse_pts(1,:), ellipse_pts(2,:), 'r-', 'LineWidth', 2);

xlabel('Dimension 2');
ylabel('Dimension 3');
title(sprintf('Vertex %d: Dim 2 vs 3', vertex_id));
legend('Estimates (100 reps)', 'True position', 'Empirical mean', ...
    '95% ellipse (plugin)', 'Location', 'best');
grid on;
hold off;

% Plots 4-6: Marginal histograms with theoretical densities
for dim = 1:3
    subplot(2, 3, 3 + dim);
    hold on;

    % Histogram
    histogram(X_estimates(:, dim), 15, 'Normalization', 'pdf', ...
        'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    % Theoretical normal density
    mu_theoretical = X_true(dim);
    sigma_theoretical = sqrt(Cov_theoretical(dim, dim));
    x_grid = linspace(mu_theoretical - 4*sigma_theoretical, ...
                      mu_theoretical + 4*sigma_theoretical, 200);
    y_theoretical = normpdf(x_grid, mu_theoretical, sigma_theoretical);
    plot(x_grid, y_theoretical, 'r-', 'LineWidth', 2);

    % Empirical normal density (for comparison)
    mu_empirical = X_mean(dim);
    sigma_empirical = sqrt(X_cov_empirical(dim, dim));
    y_empirical = normpdf(x_grid, mu_empirical, sigma_empirical);
    plot(x_grid, y_empirical, 'g--', 'LineWidth', 2);

    % Mark true value
    yl = ylim;
    plot([X_true(dim), X_true(dim)], yl, 'r--', 'LineWidth', 1.5);

    xlabel(sprintf('Dimension %d', dim));
    ylabel('Density');
    title(sprintf('Marginal Distribution: Dim %d', dim));
    legend('Empirical (100 reps)', 'Theoretical (plugin)', 'Empirical fit', ...
        'True value', 'Location', 'best');
    grid on;
    hold off;
end

sgtitle(sprintf('cgrdpg Vertex %d: Empirical vs Theoretical Distribution (Plugin Coverage: %.1f%%)', ...
    vertex_id, 100 * empirical_coverage), 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
output_dir = 'vertex_distribution_plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_file = fullfile(output_dir, sprintf('vertex_%d_distribution.png', vertex_id));
saveas(gcf, output_file);
fprintf('\nFigure saved to: %s\n', output_file);

%% Additional diagnostic: Q-Q plots
figure('Position', [100, 100, 1200, 400]);
for dim = 1:3
    subplot(1, 3, dim);

    % Standardize the estimates
    estimates_standardized = (X_estimates(:, dim) - X_true(dim)) / sqrt(Cov_theoretical(dim, dim));

    % Q-Q plot
    qqplot(estimates_standardized);
    xlabel('Theoretical Quantiles (N(0,1))');
    ylabel('Sample Quantiles');
    title(sprintf('Q-Q Plot: Dimension %d', dim));
    grid on;
end
sgtitle(sprintf('cgrdpg Vertex %d: Q-Q Plots (Testing Normality)', vertex_id), ...
    'FontSize', 14, 'FontWeight', 'bold');

output_file_qq = fullfile(output_dir, sprintf('vertex_%d_qq_plots.png', vertex_id));
saveas(gcf, output_file_qq);
fprintf('Q-Q plots saved to: %s\n', output_file_qq);

fprintf('\nDiagnosis:\n');
fprintf('  If empirical distribution is WIDER than theoretical ellipse:\n');
fprintf('    → Fisher information underestimates variance (too optimistic)\n');
fprintf('    → Results in lower coverage than 95%%\n\n');
fprintf('  If empirical distribution is NARROWER than theoretical ellipse:\n');
fprintf('    → Fisher information overestimates variance (too conservative)\n');
fprintf('    → Results in higher coverage than 95%%\n\n');
fprintf('  Current plugin coverage for this vertex: %.1f%%\n', 100 * empirical_coverage);
