% Run 100 replications of 3D CGRDPG simulation with fminunc - PARALLEL VERSION
% Uses MATLAB parfor for parallel execution on HPC

clear; clc;

%% Parameters (matching R simulations)
n = 1000;  % number of nodes
p_cov = 1000;  % number of covariates
d = 3;  % embedding dimension
p = 2;  % positive signatures (q = 1 negative)
tau = 0.001;
n_reps = 100;

fprintf('========================================\n');
fprintf('100 Replications: 3D CGRDPG with fminunc (PARALLEL)\n');
fprintf('========================================\n');
fprintf('n = %d, p_cov = %d, d = %d, tau = %.6f\n', n, p_cov, d, tau);
fprintf('Number of replications: %d\n\n', n_reps);

%% Setup parallel pool
% Get number of workers from environment or use default
n_workers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(n_workers) || n_workers < 1
    n_workers = min(maxNumCompThreads, 48);  % Use available cores, max 48
end

fprintf('Setting up parallel pool with %d workers...\n', n_workers);
pool = parpool('local', n_workers);
fprintf('Parallel pool ready!\n\n');

%% Configure fminunc options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display', 'off', ...  % Suppress per-iteration output
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', true);

%% Storage for results (can't use regular arrays in parfor)
results_all = cell(n_reps, 1);

%% Signature matrix (shared)
S = diag([1, 1, -1]);

%% Run replications in parallel
fprintf('Starting parallel replications...\n');
fprintf('Start time: %s\n\n', datestr(now));

tic;
parfor rep = 1:n_reps
    % Set random seed for reproducibility
    rng_state = RandStream('mt19937ar', 'Seed', rep + 1000);

    try
        %% Generate synthetic data
        % True latent positions (3D: 2 positive, 1 negative)
        X_true = randn(rng_state, n, d) * 0.5;

        % True covariate coefficients
        Z_true = randn(rng_state, p_cov, d) * 0.3;

        % Generate probabilities
        Y_true = X_true * S;
        P_net = Y_true * Y_true';
        P_net = max(min(P_net, 0.99), 0.01);  % Clamp to valid range

        % Generate adjacency matrix
        A = double(rand(rng_state, n) < P_net);
        A = triu(A, 1);  % Upper triangle only (undirected)
        A = A + A';  % Make symmetric

        % Generate covariates
        B = Z_true * X_true' + randn(rng_state, p_cov, n) * 0.1;

        %% Compute ASE baseline
        A_aug = A;
        deg = sum(A, 2);
        A_aug(1:n+1:end) = deg / (n - 1);

        [V, D] = eig(A_aug);
        [eigvals, idx] = sort(diag(D), 'descend');
        V = V(:, idx);
        X_ase = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

        % Procrustes alignment for ASE
        [~, X_ase_aligned] = procrustes(X_true, X_ase);
        sse_ase = sum((X_ase_aligned(:) - X_true(:)).^2);

        %% Fit CGRDPG with fminunc
        rep_start_time = tic;
        [X_opt, Z_opt, fval, exitflag, output] = fit_grdpg_fminunc(A, B, d, p, tau, options);
        rep_time = toc(rep_start_time);

        % Procrustes alignment for fminunc
        [~, X_opt_aligned] = procrustes(X_true, X_opt);
        sse_fminunc = sum((X_opt_aligned(:) - X_true(:)).^2);

        % Store results
        results_all{rep} = struct(...
            'rep', rep, ...
            'sse_ase', sse_ase, ...
            'sse_fminunc', sse_fminunc, ...
            'converged', (exitflag > 0), ...
            'iterations', output.iterations, ...
            'final_objective', fval, ...
            'exitflag', exitflag, ...
            'time', rep_time);

        % Print progress (approximate, may be out of order due to parallelism)
        if mod(rep, 10) == 0
            fprintf('  Completed rep %d\n', rep);
        end

    catch ME
        % Handle errors gracefully
        results_all{rep} = struct(...
            'rep', rep, ...
            'error', ME.message, ...
            'sse_ase', NaN, ...
            'sse_fminunc', NaN);
    end
end

total_time = toc;
fprintf('\nAll replications completed!\n');
fprintf('Total time: %.2f minutes (%.2f hours)\n', total_time/60, total_time/3600);
fprintf('Average time per rep: %.2f seconds\n\n', total_time/n_reps);

% Shutdown parallel pool
delete(pool);

%% Extract results from cell array
sse_ase = zeros(n_reps, 1);
sse_fminunc = zeros(n_reps, 1);
converged = zeros(n_reps, 1);
iterations = zeros(n_reps, 1);
time_elapsed = zeros(n_reps, 1);
final_objective = zeros(n_reps, 1);

for rep = 1:n_reps
    if isfield(results_all{rep}, 'error')
        sse_ase(rep) = NaN;
        sse_fminunc(rep) = NaN;
        converged(rep) = 0;
        iterations(rep) = 0;
        time_elapsed(rep) = 0;
        final_objective(rep) = NaN;
    else
        sse_ase(rep) = results_all{rep}.sse_ase;
        sse_fminunc(rep) = results_all{rep}.sse_fminunc;
        converged(rep) = results_all{rep}.converged;
        iterations(rep) = results_all{rep}.iterations;
        time_elapsed(rep) = results_all{rep}.time;
        final_objective(rep) = results_all{rep}.final_objective;
    end
end

%% Summary statistics
fprintf('========================================\n');
fprintf('Summary Statistics\n');
fprintf('========================================\n\n');

fprintf('ASE:\n');
fprintf('  Mean SSE:   %.4f\n', mean(sse_ase));
fprintf('  Median SSE: %.4f\n', median(sse_ase));
fprintf('  Std SSE:    %.4f\n', std(sse_ase));
fprintf('  Min SSE:    %.4f\n', min(sse_ase));
fprintf('  Max SSE:    %.4f\n\n', max(sse_ase));

valid_fminunc = ~isnan(sse_fminunc);
fprintf('fminunc:\n');
fprintf('  Successful runs: %d/%d (%.1f%%)\n', sum(valid_fminunc), n_reps, 100*mean(valid_fminunc));
fprintf('  Mean SSE:   %.4f\n', mean(sse_fminunc(valid_fminunc)));
fprintf('  Median SSE: %.4f\n', median(sse_fminunc(valid_fminunc)));
fprintf('  Std SSE:    %.4f\n', std(sse_fminunc(valid_fminunc)));
fprintf('  Min SSE:    %.4f\n', min(sse_fminunc(valid_fminunc)));
fprintf('  Max SSE:    %.4f\n\n', max(sse_fminunc(valid_fminunc)));

fprintf('Convergence:\n');
fprintf('  Converged: %d/%d (%.1f%%)\n', sum(converged), n_reps, 100*mean(converged));
fprintf('  Mean iterations: %.1f\n', mean(iterations(converged > 0)));
fprintf('  Mean time per rep: %.2f seconds\n\n', mean(time_elapsed(valid_fminunc)));

%% Save results
summary = struct();
summary.params = struct('n', n, 'p_cov', p_cov, 'd', d, 'p', p, 'tau', tau, 'n_reps', n_reps);
summary.sse_ase = sse_ase;
summary.sse_fminunc = sse_fminunc;
summary.converged = converged;
summary.iterations = iterations;
summary.time_elapsed = time_elapsed;
summary.final_objective = final_objective;
summary.total_time = total_time;
summary.n_workers = n_workers;

% Save summary
save('results/fminunc_100reps_n1000_parallel_summary.mat', 'summary');
fprintf('Summary saved to: results/fminunc_100reps_n1000_parallel_summary.mat\n');

% Save detailed results
save('results/fminunc_100reps_n1000_parallel_detailed.mat', 'results_all', '-v7.3');
fprintf('Detailed results saved to: results/fminunc_100reps_n1000_parallel_detailed.mat\n');

% Save individual replications
if ~exist('results/replications_parallel', 'dir')
    mkdir('results/replications_parallel');
end

fprintf('Saving individual replications...');
for rep = 1:n_reps
    rep_data = results_all{rep};
    save(sprintf('results/replications_parallel/rep_%03d.mat', rep), 'rep_data');
end
fprintf(' Done!\n\n');

fprintf('========================================\n');
fprintf('All results saved!\n');
fprintf('End time: %s\n', datestr(now));
fprintf('========================================\n');
