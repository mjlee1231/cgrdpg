function visualize_single_vertex_ase(vertex_id, results_dir, n_reps)
% VISUALIZE_SINGLE_VERTEX_ASE Create diagnostic plots for ASE estimates
%
% Inputs:
%   vertex_id   - Vertex index (1-1000)
%   results_dir - Directory containing result .mat files
%   n_reps      - Number of replications (default: 100)

if nargin < 3
    n_reps = 100;
end

d = 3;
alpha = 0.05;
chi2_crit = chi2inv(1 - alpha, d);
eps_clip = 1e-10;

% Storage for estimates across reps
X_estimates = nan(n_reps, d);
X_true = nan(1, d);
prec_ase_plugin = nan(d, d);
coverage_plugin = nan(n_reps, 1);

% Load all replications
n_loaded = 0;
for rep = 1:n_reps
    filename = fullfile(results_dir, sprintf('rep_%03d.mat', rep));

    if exist(filename, 'file')
        data = load(filename);
        X_estimates(rep, :) = data.X_ase(vertex_id, :);

        if rep == 1
            X_true = data.X0(vertex_id, :);
            X_ase_plugin = data.X_ase;
            S_estimated = data.S_estimated;
            % Compute plugin precision matrix for this vertex
            prec_ase_plugin = compute_prec_ase(vertex_id, X_ase_plugin, S_estimated, eps_clip);
        end

        % Check coverage for this vertex using plugin precision
        Prec_plugin_rep = compute_prec_ase(vertex_id, data.X_ase, data.S_estimated, eps_clip);
        err = data.X_ase(vertex_id, :) - X_true;
        dist = err * Prec_plugin_rep * err';  % Note: no (n+p_cov) scaling for ASE
        coverage_plugin(rep) = (dist <= chi2_crit);

        n_loaded = n_loaded + 1;
    end
end

if n_loaded == 0
    warning('No replications found for vertex %d', vertex_id);
    return;
end

% Compute empirical statistics
X_mean = mean(X_estimates, 1, 'omitnan');
X_cov_empirical = cov(X_estimates, 'omitrows');
empirical_coverage = mean(coverage_plugin, 'omitnan');

% Theoretical covariance from plugin precision matrix
% For ASE: Cov = Prec^{-1} (no scaling factor)
Cov_theoretical = inv(prec_ase_plugin);

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

sgtitle(sprintf('ASE Vertex %d: Empirical vs Theoretical Distribution (Plugin Coverage: %.1f%%)', ...
    vertex_id, 100 * empirical_coverage), 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
output_dir = 'vertex_distribution_plots_ase';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_file = fullfile(output_dir, sprintf('vertex_%d_ase_distribution.png', vertex_id));
saveas(gcf, output_file);

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
sgtitle(sprintf('ASE Vertex %d: Q-Q Plots (Testing Normality)', vertex_id), ...
    'FontSize', 14, 'FontWeight', 'bold');

output_file_qq = fullfile(output_dir, sprintf('vertex_%d_ase_qq_plots.png', vertex_id));
saveas(gcf, output_file_qq);

close all;

end
