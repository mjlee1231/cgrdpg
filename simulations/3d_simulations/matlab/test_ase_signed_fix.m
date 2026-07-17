% Quick test: verify X_ase_signed matches R after fix
clear; clc;
addpath('core');

fprintf('Testing ASE signed fix\n');
fprintf('======================\n\n');

% Parameters
n = 1000;
d = 3;
p = 2;

% Set same seed as R test
rng(599);

% Generate same random network as R test
A = double(rand(n) < 0.3);
A = triu(A, 1);
A = A + A';

% Augmented adjacency
deg = sum(A, 2);
A(1:n+1:end) = deg / (n - 1);

% Fit ASE
[X_ase_unsigned, X_ase_signed, S_estimated] = fit_ase(A, d, p);

fprintf('S_estimated = diag([%d, %d, %d])\n\n', diag(S_estimated));

% Check if X_signed == X * S
X_signed_check = X_ase_unsigned * S_estimated;
diff = max(abs(X_ase_signed(:) - X_signed_check(:)));

fprintf('Max difference between X_signed and X*S: %.2e\n', diff);

if diff < 1e-10
    fprintf('✓ PASS: X_signed == X * S_estimated\n\n');
else
    fprintf('✗ FAIL: X_signed != X * S_estimated\n\n');
end

% Compare first row with R output
fprintf('First row comparison with R:\n');
fprintf('  X[1,:] (unsigned):  [%.9f, %.9f, %.9f]\n', X_ase_unsigned(1,:));
fprintf('  X_signed[1,:]:      [%.9f, %.9f, %.9f]\n', X_ase_signed(1,:));
fprintf('\n  R expected X[1,:]:  [-0.573896819,  0.076413091, -0.001289974]\n');
fprintf('  R expected X_sgn:   [-0.573896819,  0.076413091,  0.001289974]\n');
fprintf('\nNote: Signs may differ due to eigenvector ambiguity, but pattern should match.\n');
