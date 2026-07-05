% Compare MATLAB fminunc results vs R Fisher scoring results
% Analyzes SSE distributions and convergence properties

clear; clc;

fprintf('========================================\n');
fprintf('Comparing MATLAB fminunc vs R results\n');
fprintf('========================================\n\n');

%% Load MATLAB results
if ~exist('results/fminunc_100reps_n1000_summary.mat', 'file')
    error('MATLAB results not found. Run run_100_replications.m first.');
end

matlab_data = load('results/fminunc_100reps_n1000_summary.mat');
summary = matlab_data.summary;

fprintf('MATLAB Results Loaded:\n');
fprintf('  Replications: %d\n', summary.params.n_reps);
fprintf('  n = %d, d = %d\n\n', summary.params.n, summary.params.d);

%% Display MATLAB statistics
fprintf('MATLAB fminunc Statistics:\n');
fprintf('----------------------------------------\n');
fprintf('ASE:\n');
fprintf('  Mean SSE:   %.4f\n', mean(summary.sse_ase));
fprintf('  Median SSE: %.4f\n', median(summary.sse_ase));
fprintf('  Std SSE:    %.4f\n', std(summary.sse_ase));
fprintf('  Range:      [%.4f, %.4f]\n\n', min(summary.sse_ase), max(summary.sse_ase));

valid_idx = ~isnan(summary.sse_fminunc);
fprintf('fminunc (MSLE):\n');
fprintf('  Successful: %d/%d (%.1f%%)\n', sum(valid_idx), length(valid_idx), 100*mean(valid_idx));
fprintf('  Mean SSE:   %.4f\n', mean(summary.sse_fminunc(valid_idx)));
fprintf('  Median SSE: %.4f\n', median(summary.sse_fminunc(valid_idx)));
fprintf('  Std SSE:    %.4f\n', std(summary.sse_fminunc(valid_idx)));
fprintf('  Range:      [%.4f, %.4f]\n\n', min(summary.sse_fminunc(valid_idx)), max(summary.sse_fminunc(valid_idx)));

fprintf('Convergence:\n');
fprintf('  Converged:  %d/%d (%.1f%%)\n', sum(summary.converged), length(summary.converged), 100*mean(summary.converged));
fprintf('  Mean iterations: %.1f\n', mean(summary.iterations(summary.converged > 0)));
fprintf('  Mean time: %.2f sec\n\n', mean(summary.time_elapsed(valid_idx)));

%% Try to load and compare with R results
% Note: User needs to create R results comparison files
% This section provides template for comparison

fprintf('========================================\n');
fprintf('R Results Comparison (if available)\n');
fprintf('========================================\n\n');

fprintf('To compare with R results:\n');
fprintf('1. From R simulation results, extract SSE values:\n');
fprintf('   R> library(R.matlab)\n');
fprintf('   R> sse_r <- sapply(results_list, function(x) x$sse)\n');
fprintf('   R> writeMat("r_sse_100reps.mat", sse_ase=sse_r["ase",],\n');
fprintf('               sse_ose=sse_r["ose",], sse_cgrdpg=sse_r["cgrdpg",])\n\n');

if exist('results/r_sse_100reps.mat', 'file')
    fprintf('R results found! Loading...\n');
    r_data = load('results/r_sse_100reps.mat');

    fprintf('\nR Fisher Scoring Statistics:\n');
    fprintf('----------------------------------------\n');

    fprintf('ASE (R):\n');
    fprintf('  Mean SSE:   %.4f\n', mean(r_data.sse_ase));
    fprintf('  Median SSE: %.4f\n', median(r_data.sse_ase));
    fprintf('  Std SSE:    %.4f\n', std(r_data.sse_ase));
    fprintf('  Range:      [%.4f, %.4f]\n\n', min(r_data.sse_ase), max(r_data.sse_ase));

    fprintf('OSE (R):\n');
    fprintf('  Mean SSE:   %.4f\n', mean(r_data.sse_ose));
    fprintf('  Median SSE: %.4f\n', median(r_data.sse_ose));
    fprintf('  Std SSE:    %.4f\n', std(r_data.sse_ose));
    fprintf('  Range:      [%.4f, %.4f]\n\n', min(r_data.sse_ose), max(r_data.sse_ose));

    fprintf('CGRDPG/Fisher (R):\n');
    fprintf('  Mean SSE:   %.4f\n', mean(r_data.sse_cgrdpg));
    fprintf('  Median SSE: %.4f\n', median(r_data.sse_cgrdpg));
    fprintf('  Std SSE:    %.4f\n', std(r_data.sse_cgrdpg));
    fprintf('  Range:      [%.4f, %.4f]\n\n', min(r_data.sse_cgrdpg), max(r_data.sse_cgrdpg));

    % Create comparison plots
    figure('Position', [100, 100, 1200, 400]);

    % Histogram comparison
    subplot(1, 3, 1);
    hold on;
    histogram(summary.sse_ase, 20, 'FaceAlpha', 0.5, 'DisplayName', 'ASE (MATLAB)');
    histogram(r_data.sse_ase, 20, 'FaceAlpha', 0.5, 'DisplayName', 'ASE (R)');
    xlabel('SSE'); ylabel('Frequency');
    title('ASE Comparison');
    legend('Location', 'best');
    grid on;

    subplot(1, 3, 2);
    hold on;
    histogram(r_data.sse_ose, 20, 'FaceAlpha', 0.5, 'DisplayName', 'OSE (R)');
    xlabel('SSE'); ylabel('Frequency');
    title('OSE (R only - unstable)');
    legend('Location', 'best');
    grid on;

    subplot(1, 3, 3);
    hold on;
    histogram(summary.sse_fminunc(valid_idx), 20, 'FaceAlpha', 0.5, 'DisplayName', 'fminunc (MATLAB)');
    histogram(r_data.sse_cgrdpg, 20, 'FaceAlpha', 0.5, 'DisplayName', 'Fisher (R)');
    xlabel('SSE'); ylabel('Frequency');
    title('MSLE: MATLAB vs R');
    legend('Location', 'best');
    grid on;

    saveas(gcf, 'results/matlab_vs_r_comparison.png');
    fprintf('Comparison plot saved to: results/matlab_vs_r_comparison.png\n\n');

    % Statistical comparison
    fprintf('========================================\n');
    fprintf('Direct Comparison\n');
    fprintf('========================================\n\n');

    fprintf('MSLE Methods:\n');
    fprintf('  MATLAB fminunc mean SSE: %.4f\n', mean(summary.sse_fminunc(valid_idx)));
    fprintf('  R Fisher mean SSE:       %.4f\n', mean(r_data.sse_cgrdpg));
    fprintf('  Improvement:             %.2f%%\n\n', ...
        100 * (mean(r_data.sse_cgrdpg) - mean(summary.sse_fminunc(valid_idx))) / mean(r_data.sse_cgrdpg));

    fprintf('Stability (Max SSE):\n');
    fprintf('  MATLAB fminunc: %.4f\n', max(summary.sse_fminunc(valid_idx)));
    fprintf('  R Fisher:       %.4f\n', max(r_data.sse_cgrdpg));
    fprintf('  R OSE:          %.4f (UNSTABLE!)\n\n', max(r_data.sse_ose));

else
    fprintf('R results not found. Skipping comparison.\n');
    fprintf('See instructions above to create r_sse_100reps.mat\n\n');
end

%% Create MATLAB-only distribution plot
figure('Position', [100, 100, 800, 600]);
hold on;
histogram(summary.sse_ase, 30, 'FaceAlpha', 0.6, 'DisplayName', 'ASE');
histogram(summary.sse_fminunc(valid_idx), 30, 'FaceAlpha', 0.6, 'DisplayName', 'fminunc (MSLE)');
xlabel('SSE', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('SSE Distribution: MATLAB fminunc (n=1000, d=3, 100 reps)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);

saveas(gcf, 'results/matlab_sse_distribution.png');
fprintf('MATLAB distribution plot saved to: results/matlab_sse_distribution.png\n');

fprintf('\n========================================\n');
fprintf('Analysis complete!\n');
fprintf('========================================\n');
