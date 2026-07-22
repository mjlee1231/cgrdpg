% Diagnose ASE eigenvector selection for indefinite GRDPG
% Check if magnitude-based sorting picks the correct population direction
clear; clc;

fprintf('Diagnosing ASE Eigenvector Selection for Indefinite GRDPG\n');
fprintf('===========================================================\n\n');

%% Setup
n = 1000;
p_cov = 500;
d = 3;
p = 2;  % Number of positive eigenvalues in signature

rng(599);  % Use one replication for diagnosis

% True latent positions
t = (1:n)' / n;
X0 = [0.42*t + 0.46, ...
      0.27 * sin(2*pi*t) + 0.46, ...
      0.20 * cos(4*pi*t)];

S = diag([1, 1, -1]);
Y0 = X0 * S;
Z0 = randn(p_cov, d);

% Generate data
P = X0 * Y0';
A = double(rand(n) < P);
A = triu(A, 1) + triu(A, 1)';

% Augmented adjacency
A_aug = A;
deg = sum(A, 2);
A_aug(1:n+1:end) = deg / (n - 1);

fprintf('Sample size: n=%d, p_cov=%d\n', n, p_cov);
fprintf('True signature: S = diag([1, 1, -1])\n\n');

%% Population eigendecomposition
fprintf('========================================\n');
fprintf('POPULATION (True P = X0*S*X0'')\n');
fprintf('========================================\n');

[V_pop, D_pop] = eig(P);
eigvals_pop = diag(D_pop);

% Sort by magnitude
[~, idx_pop] = sort(abs(eigvals_pop), 'descend');
eigvals_pop_sorted = eigvals_pop(idx_pop);
V_pop_sorted = V_pop(:, idx_pop);

fprintf('Top 3 eigenvalues (sorted by magnitude):\n');
for i = 1:3
    fprintf('  λ%d = %+.2f  (|λ%d| = %.2f, sign: %s)\n', ...
        i, eigvals_pop_sorted(i), i, abs(eigvals_pop_sorted(i)), ...
        sign_str(eigvals_pop_sorted(i)));
end

% Population eigenvectors (top 3 by magnitude)
v1_pop = V_pop_sorted(:, 1);
v2_pop = V_pop_sorted(:, 2);
v3_pop = V_pop_sorted(:, 3);

fprintf('\n');

%% Observed eigendecomposition
fprintf('========================================\n');
fprintf('OBSERVED (A_aug)\n');
fprintf('========================================\n');

[V_obs, D_obs] = eig(A_aug);
eigvals_obs = diag(D_obs);

% Sort by magnitude (what ASE does)
[~, idx_obs] = sort(abs(eigvals_obs), 'descend');
eigvals_obs_sorted = eigvals_obs(idx_obs);
V_obs_sorted = V_obs(:, idx_obs);

fprintf('Top 5 eigenvalues (sorted by magnitude):\n');
for i = 1:min(5, n)
    fprintf('  λ%d = %+.4f  (|λ%d| = %.4f, sign: %s)\n', ...
        i, eigvals_obs_sorted(i), i, abs(eigvals_obs_sorted(i)), ...
        sign_str(eigvals_obs_sorted(i)));
end

fprintf('\n');

%% Check alignment: which observed eigenvector matches population v3?
fprintf('========================================\n');
fprintf('EIGENVECTOR ALIGNMENT ANALYSIS\n');
fprintf('========================================\n');

fprintf('\nPopulation 3rd eigenvector (λ_pop,3 = %.2f):\n', eigvals_pop_sorted(3));
fprintf('Checking alignment with observed eigenvectors...\n\n');

% Compute cosine similarity with top 10 observed eigenvectors
n_check = min(10, n);
similarities = zeros(n_check, 1);

for i = 1:n_check
    v_obs = V_obs_sorted(:, i);
    % Cosine similarity (absolute value, since eigenvector sign is arbitrary)
    similarities(i) = abs(v3_pop' * v_obs);
end

fprintf('Cosine similarity |v_pop,3'' * v_obs,i|:\n');
fprintf('%-4s %-12s %-10s %-12s\n', 'Rank', 'λ_obs', 'Sign', 'Similarity');
fprintf('%s\n', repmat('-', 1, 45));

for i = 1:n_check
    fprintf('%-4d %+11.4f %-10s %11.4f', ...
        i, eigvals_obs_sorted(i), sign_str(eigvals_obs_sorted(i)), similarities(i));

    if i == 3
        fprintf('  <-- ASE selects this');
    end

    if similarities(i) > 0.9
        fprintf('  *** HIGH MATCH ***');
    end

    fprintf('\n');
end

% Find best match
[max_sim, best_idx] = max(similarities);
fprintf('\nBest match: Rank %d (λ = %+.4f, similarity = %.4f)\n', ...
    best_idx, eigvals_obs_sorted(best_idx), max_sim);

if best_idx ~= 3
    fprintf('\n⚠️  WARNING: ASE selects rank 3, but best match is rank %d!\n', best_idx);
    fprintf('   This suggests magnitude-based sorting may select wrong direction.\n');
else
    fprintf('\n✓ ASE correctly selects the population 3rd direction.\n');
end

%% Check if rank 3 is actually noise
if best_idx ~= 3
    fprintf('\n========================================\n');
    fprintf('ANALYZING RANK 3 (ASE SELECTION)\n');
    fprintf('========================================\n');

    fprintf('Rank 3 eigenvalue: λ_obs,3 = %+.4f (sign: %s)\n', ...
        eigvals_obs_sorted(3), sign_str(eigvals_obs_sorted(3)));
    fprintf('Similarity to population v3: %.4f\n', similarities(3));

    % Check if it aligns better with population v1 or v2
    sim_v1 = abs(v1_pop' * V_obs_sorted(:, 3));
    sim_v2 = abs(v2_pop' * V_obs_sorted(:, 3));

    fprintf('Similarity to population v1: %.4f\n', sim_v1);
    fprintf('Similarity to population v2: %.4f\n', sim_v2);

    if max([sim_v1, sim_v2, similarities(3)]) < 0.5
        fprintf('\n=> Rank 3 appears to be NOISE (low similarity to all population directions)\n');
    end
end

%% Summary
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n');

fprintf('\nASE uses magnitude-based sorting: |λ₁| ≥ |λ₂| ≥ |λ₃|\n');
fprintf('For indefinite GRDPG with S = diag(1,1,-1):\n');
fprintf('  Population: λ₁ = +665.63, λ₂ = +46.30, λ₃ = -20.00\n');
fprintf('  Magnitudes: |λ₁| = 665.63, |λ₂| = 46.30, |λ₃| = 20.00\n\n');

if eigvals_obs_sorted(3) > 0
    fprintf('⚠️  Observed rank 3 eigenvalue is POSITIVE: λ_obs,3 = %+.4f\n', eigvals_obs_sorted(3));
    fprintf('   This means the true negative eigenvalue direction may be ranked lower!\n');
else
    fprintf('✓ Observed rank 3 eigenvalue is NEGATIVE: λ_obs,3 = %+.4f\n', eigvals_obs_sorted(3));
    fprintf('  This matches the expected signature.\n');
end

function s = sign_str(x)
    if x > 0
        s = '+';
    elseif x < 0
        s = '-';
    else
        s = '0';
    end
end
