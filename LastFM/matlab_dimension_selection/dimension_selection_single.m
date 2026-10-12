function dimension_selection_single(d)
% DIMENSION_SELECTION_SINGLE Fit GRDPG for a single dimension value
%
% Usage: dimension_selection_single(d)
%   d - embedding dimension (1, 2, 3, 4, or 5)
%
% Outputs saved to: results_matlab/results_d{d}.mat

fprintf('============================================================================\n');
fprintf('  DIMENSION SELECTION: d = %d\n', d);
fprintf('  MATLAB Implementation\n');
fprintf('============================================================================\n\n');

fprintf('Job started: %s\n', datetime);
fprintf('Dimension: d = %d\n\n', d);

addpath('core');  % Add core functions to path

%% STEP 1: Load and preprocess data
fprintf('STEP 1: Loading and preprocessing data...\n');

% Load data
load('../CovariatesMatrix.mat');  % Load MATLAB format (created by prepare_data_for_matlab.R)
A_edges = readmatrix('../lastfm_asia_edges.csv');
labels_data = readmatrix('../lastfm_asia_target.csv');

% Extract target column (last column if CSV has multiple columns)
if size(labels_data, 2) > 1
    labels = labels_data(:, end);  % Take last column (target)
else
    labels = labels_data;
end

% Convert edgelist to adjacency matrix
n_nodes = max(A_edges(:)) + 1;  % +1 because edges are 0-indexed in CSV
A = zeros(n_nodes, n_nodes);
for i = 1:size(A_edges, 1)
    A(A_edges(i,1)+1, A_edges(i,2)+1) = 1;  % +1 for MATLAB 1-indexing
end
A = A + A';  % Symmetrize

% Debug: Check dimensions
fprintf('  Loaded data: A is %dx%d, labels is %dx%d, sparseX is %dx%d\n', ...
    size(A, 1), size(A, 2), size(labels, 1), size(labels, 2), ...
    size(sparseX, 1), size(sparseX, 2));

% Preprocessing (match R version)
label = labels + 1;  % Convert to 1-indexed
ind5 = label ~= 5;

% Debug: Check filtering
fprintf('  Filtering: keeping %d/%d nodes (excluding label 5)\n', sum(ind5), length(ind5));
A = A(ind5, ind5);
X = sparseX(ind5, :);
label = label(ind5);

% Remove zero-sum columns
X = X(:, sum(X) > 0);

% Adjust labels
labelnew = label;
label(labelnew > 5) = label(labelnew > 5) - 1;

% Select class sizes
sizes = histcounts(label, 'BinMethod', 'integers');
class_indices = find(sizes > 300 & sizes < 1000);
ind_select = ismember(label, class_indices);

Aselect = A(ind_select, ind_select);
Xselect = X(ind_select, :);
labelselect = label(ind_select);

% Artist filtering
dartist_all = sum(X);
dartist = sum(Xselect);
prop = dartist ./ dartist_all;
n_select = size(Aselect, 1);
prob = 1 - round(min(n_select/2, 600)/size(Xselect, 2), 4);
Xselect = Xselect(:, prop > quantile(prop, prob));

% Node filtering
dnode = sum(Xselect, 2);
like = dnode > 0;
Xselect = Xselect(like, :);
Aselect = Aselect(like, like);
labelselect = labelselect(like);

% Remove low-degree nodes
dselect = sum(Aselect, 2);
while sum(dselect <= 1) > 0
    keep = dselect > 1;
    Aselect = Aselect(keep, keep);
    Xselect = Xselect(keep, :);
    labelselect = labelselect(keep);
    dselect = sum(Aselect, 2);
end

n = size(Aselect, 1);
p_cov = size(Xselect, 2);
B = Xselect';  % Transpose to match R convention

fprintf('  n=%d, p=%d, k=%d countries\n', n, p_cov, length(unique(labelselect)));

% Compute tau
P_temp = mean(Aselect(triu(true(n), 1)));
tau = min(P_temp / 10, 0.001);
fprintf('  Network density: %.4f, tau: %.6f\n\n', P_temp, tau);

%% STEP 2: Fit model
fprintf('============================================================================\n');
fprintf('STEP 2: Fitting model with d = %d\n', d);
fprintf('============================================================================\n\n');

% Set up optimization options
options = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective');

fprintf('Fit started: %s\n\n', datetime);
tic;

% Fit model
p_sig = d;  % All positive signature for simplicity (can adjust)
[X_opt, Z_opt, fval, exitflag, output, S_estimated] = ...
    fit_grdpg_fminunc_surrogate(Aselect, B, d, p_sig, tau, options);

elapsed = toc;

converged = (exitflag == 1);
iters = output.iterations;

fprintf('\n✓ Fitting completed in %.2f seconds (%.2f hours)\n', elapsed, elapsed/3600);
fprintf('  Converged: %s\n', string(converged));
fprintf('  Iterations: %d\n', iters);
fprintf('  Final objective: %.4f\n', -fval);

%% STEP 3: Compute information criteria
fprintf('\n============================================================================\n');
fprintf('STEP 3: Computing information criteria\n');
fprintf('============================================================================\n\n');

% Negative log-likelihood
neg_log_lik = -fval;

% Number of parameters: n*d + p*d - d^2
n_params = n * d + p_cov * d - d * d;

% BIC and AIC
BIC = 2 * neg_log_lik + n_params * log(n);
AIC = 2 * neg_log_lik + 2 * n_params;

fprintf('  Neg-log-likelihood: %.4f\n', neg_log_lik);
fprintf('  Number of parameters: %d\n', n_params);
fprintf('  BIC: %.4f\n', BIC);
fprintf('  AIC: %.4f\n\n', AIC);

%% STEP 4: Save results
fprintf('============================================================================\n');
fprintf('STEP 4: Saving results\n');
fprintf('============================================================================\n\n');

result = struct();
result.d = d;
result.neg_log_lik = neg_log_lik;
result.n_params = n_params;
result.BIC = BIC;
result.AIC = AIC;
result.converged = converged;
result.iters = iters;
result.time_sec = elapsed;
result.data_info = struct('n', n, 'p', p_cov, 'k', length(unique(labelselect)));
result.latent_positions = X_opt;
result.covariate_effects = Z_opt;
result.signature = S_estimated;
result.timestamp = datetime;

% Create output directory if needed (with settings info)
results_dir = 'results_matlab_maxit30_tol0.1';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

output_file = sprintf('%s/results_d%d.mat', results_dir, d);
save(output_file, 'result');

fprintf('✓ Results saved to: %s\n\n', output_file);

%% STEP 5: Summary
fprintf('============================================================================\n');
fprintf('JOB COMPLETE\n');
fprintf('============================================================================\n\n');

fprintf('SUMMARY:\n');
fprintf('  Dimension: d = %d\n', d);
fprintf('  Converged: %s (%d iterations)\n', string(converged), iters);
fprintf('  Time: %.2f hours\n', elapsed/3600);
fprintf('  Neg-log-lik: %.4f\n', neg_log_lik);
fprintf('  BIC: %.4f\n', BIC);
fprintf('  AIC: %.4f\n', AIC);
fprintf('  Output: %s\n', output_file);
fprintf('\nJob completed: %s\n', datetime);
fprintf('============================================================================\n');

end
