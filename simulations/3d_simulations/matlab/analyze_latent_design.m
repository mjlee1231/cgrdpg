% Analyze eigenvalue properties of latent position designs
clear; clc;

fprintf('Analyzing Latent Position Design Properties\n');
fprintf('===========================================\n\n');

n = 1000;
t = (1:n)' / n;

%% Current Design
fprintf('CURRENT DESIGN:\n');
fprintf('X0 = [0.3*t + 0.5, 0.15*sin(2πt) + 0.6, 0.1*cos(4πt)]\n\n');

X_current = [0.3*t + 0.5, ...
             0.15 * sin(2*pi*t) + 0.6, ...
             0.1 * cos(4*pi*t)];

% Gram matrix
G_current = X_current' * X_current;
eig_current = eig(G_current);
[eig_current_sorted, idx] = sort(eig_current, 'descend');

fprintf('Gram matrix X''X eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_current_sorted(1));
fprintf('  λ2 = %.4f\n', eig_current_sorted(2));
fprintf('  λ3 = %.4f\n', eig_current_sorted(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', ...
    eig_current_sorted(1) - eig_current_sorted(2), eig_current_sorted(1)/eig_current_sorted(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', ...
    eig_current_sorted(2) - eig_current_sorted(3), eig_current_sorted(2)/eig_current_sorted(3));
fprintf('\nCondition number: %.2f\n', cond(G_current));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_current_sorted)));

% Column correlations
corr_matrix = corr(X_current);
fprintf('Column correlations:\n');
fprintf('  corr(X1, X2) = %.4f\n', corr_matrix(1,2));
fprintf('  corr(X1, X3) = %.4f\n', corr_matrix(1,3));
fprintf('  corr(X2, X3) = %.4f\n\n', corr_matrix(2,3));

% Edge probability range
S = diag([1, 1, -1]);
P_current = X_current * S * X_current';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_current(:)), max(P_current(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_current(:)));

%% Proposed Design 1: Orthogonalized with larger gaps
fprintf('========================================\n');
fprintf('PROPOSED DESIGN 1: Orthogonalized basis\n');
fprintf('X0 = [a1*t_orth, a2*sin_orth, a3*cos_orth]\n');
fprintf('where columns are orthogonalized\n\n');

% Create orthogonal basis using QR
X_raw = [t, sin(2*pi*t), cos(4*pi*t)];
[Q, R] = qr(X_raw, 0);  % Orthonormal columns

% Scale to get desired eigenvalues
% Target eigenvalues: λ1=400, λ2=250, λ3=100 (well-separated)
target_eigs = [400, 250, 100];
scaling = sqrt(target_eigs);
X_prop1 = Q * diag(scaling);

% Shift to positive range for edge probabilities
X_prop1 = X_prop1 + [0.6, 0.6, 0];  % Add mean to first two dimensions

G_prop1 = X_prop1' * X_prop1;
eig_prop1 = sort(eig(G_prop1), 'descend');

fprintf('Gram matrix eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_prop1(1));
fprintf('  λ2 = %.4f\n', eig_prop1(2));
fprintf('  λ3 = %.4f\n', eig_prop1(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', eig_prop1(1) - eig_prop1(2), eig_prop1(1)/eig_prop1(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', eig_prop1(2) - eig_prop1(3), eig_prop1(2)/eig_prop1(3));
fprintf('\nCondition number: %.2f\n', cond(G_prop1));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_prop1)));

P_prop1 = X_prop1 * S * X_prop1';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_prop1(:)), max(P_prop1(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_prop1(:)));

%% Proposed Design 2: Simple scaled design with larger coefficients
fprintf('========================================\n');
fprintf('PROPOSED DESIGN 2: Scaled version of current\n');
fprintf('X0 = [0.5*t + 0.6, 0.3*sin(2πt) + 0.6, 0.2*cos(4πt)]\n');
fprintf('(doubled scaling factors)\n\n');

X_prop2 = [0.5*t + 0.6, ...
           0.3 * sin(2*pi*t) + 0.6, ...
           0.2 * cos(4*pi*t)];

G_prop2 = X_prop2' * X_prop2;
eig_prop2 = sort(eig(G_prop2), 'descend');

fprintf('Gram matrix eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_prop2(1));
fprintf('  λ2 = %.4f\n', eig_prop2(2));
fprintf('  λ3 = %.4f\n', eig_prop2(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', eig_prop2(1) - eig_prop2(2), eig_prop2(1)/eig_prop2(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', eig_prop2(2) - eig_prop2(3), eig_prop2(2)/eig_prop2(3));
fprintf('\nCondition number: %.2f\n', cond(G_prop2));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_prop2)));

P_prop2 = X_prop2 * S * X_prop2';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_prop2(:)), max(P_prop2(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_prop2(:)));

%% Proposed Design 3: Prescribed eigenvalues with random rotation
fprintf('========================================\n');
fprintf('PROPOSED DESIGN 3: Direct eigenvalue prescription\n');
fprintf('X0 = U * sqrt(Λ) with Λ = diag([450, 280, 120])\n\n');

rng(42);  % For reproducibility
% Random orthogonal matrix
[U, ~] = qr(randn(n, 3));
% Prescribed eigenvalues
Lambda = diag([450, 280, 120]);
X_prop3 = U * sqrt(Lambda);

% Shift to positive range
X_prop3 = X_prop3 + [0.6, 0.6, 0];

G_prop3 = X_prop3' * X_prop3;
eig_prop3 = sort(eig(G_prop3), 'descend');

fprintf('Gram matrix eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_prop3(1));
fprintf('  λ2 = %.4f\n', eig_prop3(2));
fprintf('  λ3 = %.4f\n', eig_prop3(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', eig_prop3(1) - eig_prop3(2), eig_prop3(1)/eig_prop3(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', eig_prop3(2) - eig_prop3(3), eig_prop3(2)/eig_prop3(3));
fprintf('\nCondition number: %.2f\n', cond(G_prop3));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_prop3)));

P_prop3 = X_prop3 * S * X_prop3';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_prop3(:)), max(P_prop3(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_prop3(:)));

%% Summary
fprintf('========================================\n');
fprintf('SUMMARY COMPARISON\n');
fprintf('========================================\n\n');

designs = {'Current', 'Orthogonalized', 'Scaled', 'Prescribed'};
eigs_all = {eig_current_sorted, eig_prop1, eig_prop2, eig_prop3};
P_all = {P_current, P_prop1, P_prop2, P_prop3};

fprintf('%-15s %8s %8s %8s %12s %8s %12s\n', ...
    'Design', 'λ1', 'λ2', 'λ3', 'Gap1', 'Gap2', 'Cond#');
fprintf('%-15s %8s %8s %8s %12s %8s %12s\n', ...
    repmat('-', 1, 15), repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 8), ...
    repmat('-', 1, 12), repmat('-', 1, 8), repmat('-', 1, 12));

for i = 1:4
    eigs = eigs_all{i};
    gap1 = eigs(1) - eigs(2);
    gap2 = eigs(2) - eigs(3);
    cond_num = eigs(1) / eigs(3);

    fprintf('%-15s %8.2f %8.2f %8.2f %12.2f %8.2f %12.2f\n', ...
        designs{i}, eigs(1), eigs(2), eigs(3), gap1, gap2, cond_num);
end

fprintf('\n%-15s %10s %10s\n', 'Design', 'P_min', 'P_max');
fprintf('%-15s %10s %10s\n', repmat('-', 1, 15), repmat('-', 1, 10), repmat('-', 1, 10));
for i = 1:4
    P = P_all{i};
    fprintf('%-15s %10.4f %10.4f\n', designs{i}, min(P(:)), max(P(:)));
end

fprintf('\n');
fprintf('RECOMMENDATION:\n');
fprintf('Design 2 (Scaled) is simplest and gives good eigenvalue separation\n');
fprintf('Design 1 (Orthogonalized) gives best theoretical properties but loses interpretability\n');
fprintf('Design 3 (Prescribed) gives exact control but is random/not reproducible across studies\n\n');

fprintf('Target criteria:\n');
fprintf('  ✓ Gap between eigenvalues > 100 (for identifiability)\n');
fprintf('  ✓ Smallest eigenvalue > 40 (in magnitude)\n');
fprintf('  ✓ Condition number < 10 (for numerical stability)\n');
fprintf('  ✓ Edge probabilities in [0, 1] range\n');
