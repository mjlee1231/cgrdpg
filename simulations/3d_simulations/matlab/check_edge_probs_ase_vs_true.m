% Compare edge probabilities at TRUE positions vs ASE initialization
clear; clc;
addpath('core');

fprintf('Edge Probabilities: TRUE vs ASE initialization\n');
fprintf('================================================\n\n');

n = 1000;
p_cov = 500;
d = 3;
p = 2;
tau = 0.001;
eps_clip = 1e-10;

% Set seed
rng(598 + 1);

% Generate data (NEW design)
t = (1:n)' / n;
X0 = [0.3*t + 0.5, 0.15 * sin(2*pi*t) + 0.6, 0.1 * cos(4*pi*t)];
S = diag([1, 1, -1]);
Y0 = X0 * S;
P = X0 * Y0';

% Generate network
A = double(rand(n) < P);
A = triu(A, 1) + triu(A, 1)';
B = randn(p_cov, d) * X0' + randn(p_cov, n);

% ASE initialization
[X_ase_unsigned, X_ase_signed, S_estimated] = fit_ase(A, d, p);

fprintf('1. TRUE POSITIONS (X0)\n');
fprintf('   Position ranges:\n');
fprintf('     Dim 1: [%.4f, %.4f]\n', min(X0(:,1)), max(X0(:,1)));
fprintf('     Dim 2: [%.4f, %.4f]\n', min(X0(:,2)), max(X0(:,2)));
fprintf('     Dim 3: [%.4f, %.4f]\n\n', min(X0(:,3)), max(X0(:,3)));

% Edge probs at true positions
LinearPred_true = X0 * S * X0';
triu_idx = triu(true(n), 1);
s_true = LinearPred_true(triu_idx);
p_true = max(min(1 ./ (1 + exp(-s_true / tau)), 1 - eps_clip), eps_clip);

fprintf('   Linear predictor s_ij:\n');
fprintf('     Range: [%.4f, %.4f]\n', min(s_true), max(s_true));
fprintf('   Edge probabilities:\n');
fprintf('     Range: [%.10f, %.10f]\n', min(p_true), max(p_true));
fprintf('     # at ceiling: %d / %d (%.1f%%)\n\n', ...
    sum(p_true >= 1 - eps_clip), length(p_true), ...
    100 * sum(p_true >= 1 - eps_clip) / length(p_true));

fprintf('2. ASE INITIALIZATION (X_ase_unsigned)\n');
fprintf('   Position ranges:\n');
fprintf('     Dim 1: [%.4f, %.4f]\n', min(X_ase_unsigned(:,1)), max(X_ase_unsigned(:,1)));
fprintf('     Dim 2: [%.4f, %.4f]\n', min(X_ase_unsigned(:,2)), max(X_ase_unsigned(:,2)));
fprintf('     Dim 3: [%.4f, %.4f]\n\n', min(X_ase_unsigned(:,3)), max(X_ase_unsigned(:,3)));

% Edge probs at ASE positions
LinearPred_ase = X_ase_unsigned * S_estimated * X_ase_unsigned';
s_ase = LinearPred_ase(triu_idx);
p_ase = max(min(1 ./ (1 + exp(-s_ase / tau)), 1 - eps_clip), eps_clip);

fprintf('   Linear predictor s_ij:\n');
fprintf('     Range: [%.4f, %.4f]\n', min(s_ase), max(s_ase));
fprintf('   Edge probabilities:\n');
fprintf('     Range: [%.10f, %.10f]\n', min(p_ase), max(p_ase));
fprintf('     # at ceiling: %d / %d (%.1f%%)\n\n', ...
    sum(p_ase >= 1 - eps_clip), length(p_ase), ...
    100 * sum(p_ase >= 1 - eps_clip) / length(p_ase));

fprintf('3. COMPARISON\n');
fprintf('   At TRUE:  %.1f%% of edge probs saturated at 1.0\n', ...
    100 * sum(p_true >= 1 - eps_clip) / length(p_true));
fprintf('   At ASE:   %.1f%% of edge probs saturated at 1.0\n\n', ...
    100 * sum(p_ase >= 1 - eps_clip) / length(p_ase));

fprintf('4. DIAGNOSIS\n');
if sum(p_ase >= 1 - eps_clip) / length(p_ase) > 0.95
    fprintf('   ❌ ASE initialization has %.1f%% saturated edge probs\n', ...
        100 * sum(p_ase >= 1 - eps_clip) / length(p_ase));
    fprintf('   This causes dpsi weights to explode:\n');
    fprintf('     dpsi = 1 / (tau * p * (1-p))\n');
    fprintf('          = 1 / (%.4f * 1.0 * 1e-10)\n', tau);
    fprintf('          = %.2e\n\n', 1 / (tau * 1.0 * 1e-10));
    fprintf('   SOLUTION: Increase tau to avoid saturation\n');
    fprintf('   Try tau = 0.01 or 0.05\n');
else
    fprintf('   ✓ ASE initialization has reasonable edge prob range\n');
    fprintf('   Current tau = %.4f should work\n', tau);
end
