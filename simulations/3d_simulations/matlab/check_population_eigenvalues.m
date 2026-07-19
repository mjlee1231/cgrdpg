% Check population eigenvalues from the true GRDPG model
% Understand spectral structure and why optimization is challenging

clear; clc;

fprintf('Population Eigenvalue Analysis\n');
fprintf('==============================\n\n');

% Parameters (same as simulation)
n = 1000;
d = 3;
p = 2;
q = 1;

% Set seed (rep 1)
rng(598 + 1);

%% Generate true latent positions
t = (1:n)' / n;
X0 = [0.15 * sin(2*pi*t) + 0.6, ...
      0.15 * cos(2*pi*t) + 0.6, ...
      0.15 * cos(4*pi*t)];

S = diag([1, 1, -1]);
Y0 = X0 * S;

%% Expected adjacency matrix
P = X0 * Y0';

fprintf('Edge probability statistics:\n');
fprintf('  Min:    %.6f\n', min(P(:)));
fprintf('  Max:    %.6f\n', max(P(:)));
fprintf('  Mean:   %.6f\n', mean(P(:)));
fprintf('  Median: %.6f\n', median(P(:)));
fprintf('  Std:    %.6f\n\n', std(P(:)));

%% Eigenvalue analysis of P
[V_P, D_P] = eig(P);
eigvals_P = diag(D_P);

% Sort by magnitude
[~, idx] = sort(abs(eigvals_P), 'descend');
eigvals_P_sorted = eigvals_P(idx);

fprintf('Top 10 eigenvalues of P (by magnitude):\n');
for i = 1:min(10, length(eigvals_P_sorted))
    fprintf('  λ_%d = %+10.4f\n', i, eigvals_P_sorted(i));
end
fprintf('\n');

% Focus on top d eigenvalues
top_d = eigvals_P_sorted(1:d);
fprintf('Top %d eigenvalues (embedding dimension):\n', d);
fprintf('  λ_1 = %+10.4f\n', top_d(1));
fprintf('  λ_2 = %+10.4f\n', top_d(2));
fprintf('  λ_3 = %+10.4f\n', top_d(3));
fprintf('\n');

fprintf('Eigenvalue ratios:\n');
fprintf('  λ_2 / λ_1 = %.4f\n', abs(top_d(2)) / abs(top_d(1)));
fprintf('  λ_3 / λ_1 = %.4f\n', abs(top_d(3)) / abs(top_d(1)));
fprintf('  λ_3 / λ_2 = %.4f\n', abs(top_d(3)) / abs(top_d(2)));
fprintf('\n');

%% Eigenvalue analysis of X0' * X0 (Gram matrix)
XtX = X0' * X0;
eigvals_XtX = eig(XtX);
[~, idx] = sort(abs(eigvals_XtX), 'descend');
eigvals_XtX_sorted = eigvals_XtX(idx);

fprintf('Eigenvalues of X^T X (Gram matrix):\n');
for i = 1:d
    fprintf('  λ_%d = %+10.4f\n', i, eigvals_XtX_sorted(i));
end
fprintf('\n');

fprintf('Condition number of X^T X: %.2e\n', cond(XtX));
fprintf('Condition number of P:     %.2e\n\n', cond(P));

%% Check spectral gap
fprintf('Spectral gaps:\n');
fprintf('  Gap 1-2: %.4f\n', abs(top_d(1)) - abs(top_d(2)));
fprintf('  Gap 2-3: %.4f\n', abs(top_d(2)) - abs(top_d(3)));
fprintf('  Gap 3-4: %.4f\n', abs(top_d(3)) - abs(eigvals_P_sorted(4)));
fprintf('\n');

%% Reconstruct P from top d eigenvalues
V_d = V_P(:, idx(1:d));
Lambda_d = diag(eigvals_P_sorted(1:d));
P_approx = V_d * Lambda_d * V_d';

reconstruction_error = norm(P - P_approx, 'fro') / norm(P, 'fro');
fprintf('Rank-%d approximation error: %.6f\n\n', d, reconstruction_error);

%% Check if eigenvalue distribution explains optimization difficulty
fprintf('========================================\n');
fprintf('Interpretation for Optimization\n');
fprintf('========================================\n\n');

% Check if eigenvalues are well-separated
if abs(top_d(2)) / abs(top_d(1)) > 0.9
    fprintf('⚠ WARNING: λ_2 and λ_1 are very close (ratio %.3f)\n', ...
        abs(top_d(2)) / abs(top_d(1)));
    fprintf('  This makes spectral estimation difficult.\n\n');
end

if abs(top_d(3)) / abs(top_d(2)) > 0.5
    fprintf('⚠ WARNING: λ_3 and λ_2 are relatively close (ratio %.3f)\n', ...
        abs(top_d(3)) / abs(top_d(2)));
    fprintf('  3D embedding may have ambiguity.\n\n');
end

% Check for indefinite structure
n_positive = sum(eigvals_P_sorted > 0);
n_negative = sum(eigvals_P_sorted < 0);
fprintf('Eigenvalue signs:\n');
fprintf('  Positive: %d\n', n_positive);
fprintf('  Negative: %d\n', n_negative);
fprintf('  Zero:     %d\n\n', sum(abs(eigvals_P_sorted) < 1e-10));

if sign(top_d(1)) ~= sign(top_d(2)) || sign(top_d(2)) ~= sign(top_d(3))
    fprintf('⚠ Top 3 eigenvalues have MIXED signs:\n');
    fprintf('  This creates indefinite structure (S = diag([1,1,-1]))\n');
    fprintf('  Optimization landscape may have saddle points.\n\n');
end

%% Visualize eigenvalue spectrum
fprintf('Full eigenvalue spectrum summary:\n');
fprintf('  Total eigenvalues: %d\n', length(eigvals_P_sorted));
fprintf('  Largest magnitude: %.4f\n', max(abs(eigvals_P_sorted)));
fprintf('  Effective rank (>1%% of max): %d\n', sum(abs(eigvals_P_sorted) > 0.01 * max(abs(eigvals_P_sorted))));
fprintf('\n');

% Check for numerical issues
fprintf('Numerical stability indicators:\n');
fprintf('  Max edge prob: %.6f', max(P(:)));
if max(P(:)) > 0.99
    fprintf(' ⚠ VERY CLOSE TO 1 - causes Fisher info explosion!\n');
else
    fprintf(' (OK)\n');
end
fprintf('  Min edge prob: %.6f', min(P(:)));
if min(P(:)) < 0.01
    fprintf(' ⚠ Very close to 0\n');
else
    fprintf(' (OK)\n');
end
