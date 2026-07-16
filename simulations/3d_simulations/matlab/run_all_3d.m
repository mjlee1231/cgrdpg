% Run all 100 replications for 3D coverage
fprintf('Running 3D Coverage: 100 Replications\n');
fprintf('======================================\n\n');

addpath('core');
total_start = tic;

for rep = 1:100
    fprintf('\n--- Rep %d/100 ---\n', rep);
    try
        coverage_3d_single_rep(rep);
    catch ME
        fprintf('ERROR in rep %d: %s\n', rep, ME.message);
    end
end

fprintf('\nTotal time: %.1f min\n', toc(total_start)/60);
fprintf('Now run: matlab -batch "aggregate_coverage_3d"\n');
